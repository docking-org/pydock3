class SDI:
    @staticmethod
    def from_list(lobj):
        pass
    @staticmethod
    def from_path(pobj):
        pass
    def create_subset(self, subset_indices, subset_id=None):
        pass
    def get_subset(self, subset_id):
        pass
    # adds a new block associated with a given range of the source sdi file
    def add_block(self, bstart, bfinish):
        pass
    # returns all entries associated with the given block
    def get_block(self, idx) -> list:
        pass
    # returns a single entry from the source sdi
    def get(self, idx) -> str:

class ComputeQueue(ABC):
    @abstractmethod
    def ensure_path_is_accessible(fileobj):
        pass
    @abstractmethod
    def is_file_local(fileobj):
        pass
    @abstractmethod
    def submit(self, job, deps=None, array=None, params={})
        pass

class QBlasterBase(ABC):
    # help, of course- the most important command
    @abstractmethod
    def help(self):
        pass

    # top level submit and status for program
    @abstractmethod
    def submit(self, sdi_id=None):
        pass
    @abstractmethod
    def status(self, sdi_id=None):
        pass

    # add, subtract, and list standard input files for program
    @abstractmethod
    def add_sdi(self, sdi_path):
        pass
    @abstractmethod
    def rm_sdi(self, sdi_id) # sdi_id can be a checksum or real path to the file
        pass
    @abstractmethod
    def list_sdi(self):
        pass

    # get and set the compute queue- include checks to make sure files in base are accessible from queue
    @abstractmethod
    def set_queue(self, queue_name):
        pass
    @abstractmethod
    def get_queue(self, queue_name):
        pass

    # create an archive from the current base, deactivating it
    @abstractmethod
    def create_archive(self):
        pass
    # activate an archived base as the new active base
    @abstractmethod
    def activate_archive(self, archive_path):
        pass
    # list all archives and their properties
    @abstractmethod
    def list_archives(self):
        pass
    # destroy an archive
    @abstractmethod
    def destroy_archive(self, archive_path):
        pass
    # destroy all files in this base
    @abstractmethod
    def destroy_base(self):
        pass

class LSDBase(QBlasterBase):
    def __init__(self, cfg, filebase : ParallelJobFileBase, queue : JobQueue):
        self.filebase = filebase
        self.sdibase = filebase.dir('sdi')

    ### internal stats logic, used by jobs & status command
    def get_lsd_stats(self, sdiname):
        lsd_out_dir = self.filebase.dir(sdiname)
    def get_lsd_block_stats(self, sdiname, run):
        pass
    def get_top_stats(self):
        pass

    ### queue logic
    def set_queue(self, queue_name):
        pass
    def get_queue(self):
        pass

    ### dockfiles logic
    def set_dockfiles(self, dockfiles_path):
        pass
    def get_dockfiles(self):
        pass

    ### archive/delete logic
    def list_archives(self):
        pass
    def rm_archive(self, archivename):
        pass
    def activate_archive(self, archivename):
        pass
    def destroy_base(self):
        pass

    ### SDI logic
    def add_sdi(self, sdi_path):

        new_sdi = SDI.from_path(sdi_file)
        if self.get_sdi(new_sdi.signature):
            print("sdi already exists!")
            return False
        else:
            new_sdi.source_path = sdi_path
            sdi_dest = self.filebase.dir("sdi").file(new_sdi.signature)
            new_sdi.save(sdi_dest)
        return new_sdi.signature

    def get_sdi(self, sdi_id):
        
        sdi_file = self.sdibase.file(sdi_id)
        if sdi_file.exists():
            sdi = SDI.from_file(sdi_file)
            return sdi
        else:
            return None

    def rm_sdi(self, sdi_id):

        sdi_file = self.sdibase.file(sdi_id)
        if sdi_file.exists():
            sdi = SDI.from_file(sdi_file)
            lsd_status, lsd_stats = get_lsd_stats(sdi.signature)
            if lsd_status != "":
                print("")
            sdi.destroy()

    def list_sdi(self):
        pass

    ### submit logic
    def __submit_sdi(self, sdiname, dockfilesname):
        sdi_status = get_lsd_stats(sdiname)
        sdi_status_ready = filter(sdi_status, lambda x:x[1] == "inactive")
        sdi_ready_list = [x[0] for x in sdi_status_ready]

        run = self.filebase.new_sdi_run(sdiname, subset=sdi_ready_list)

        n_jobs = 0
        blk_id = 0
        while n_jobs < self.cfg["max_queued"]:
            self.queue.submit("lsd", params={
                "LSD_SDI" : self.filebase.get_sdi(sdiname, run=r).path,
                "LSD_DOCKFILES" : self.filebase.get_dir(dockfilesname),
            })
    def submit(self, sdiname=None, confirm=False):
        if sdiname is None:
            for sdiname in self.listsdi():
                self.__submit_sdi(sdiname)
        else:
            self.__submit_sdi(sdiname)

    ### status logic
    def __status_sdi(self)
    
    def status(self, sdiname=None):
        pass

    ###################
    # example session #
    ###################
    # qblaster lsd --base s3://mybucket/mylsdbase
    # > add-sdi z22://h17p200-h17p400
    # > add-sdi z22://h18p200-h18p400
    # > set-dockfiles ~/mydockfiles
    # > get-dockfiles
    #   +------------------------------------+----------+
    #   | path                               | checksum |
    #   +------------------------------------+----------+
    #   | s3://mybucket/mylsdbase/dockfiles  | 123fda   |
    #   +------------------------------------+----------+
    # > set-queue aws_ue1 # where are access keys set?
    # > submit
    # > set-dockfiles ~/otherdockfiles
    #   Error! Active jobs are currently using dockfiles. Wait for them to finish or cancel them.
    #
    # > set-queue slurm
    #   Error! Active jobs are currently using queue. Wait for them to finish or cancel them.
    #
    # > cancel
    #   Cancelling Jobs... Done!
    #
    # > set-dockfiles ~/otherdockfiles
    #   Warning! You are trying to change your dockfiles but some poses have already been produced from them.
    #   Choose one of the following actions:
    #   1. Archive existing work, initialize a fresh run with new dockfiles
    #   2. Destroy existing work, initialize a fresh run with new dockfiles
    #   3. Continue with existing run, just swap out dockfiles.
    #   [choose from 1, 2, or 3]: 1
    #   Archive just top poses? [y/N]: y
    #
    #   Creating archive @ s3://mybucket/mylsdbase/123fda.archive
    #
    # > set-queue slurm
    #   Warning: base in "S3" is accessible but not local to queue "slurm", expect additional data transfer fees
    #
    # > list-sdi
    #   +-----------------------+-----------+
    #   | path                  | checksum  |
    #   +-----------------------+-----------+
    #   | z22://h17p200-h17p400 | a342ve    |
    #   | z22://h18p200-h18p400 | 1usq3r    |
    #   +-----------------------+-----------+
    # > submit z22://h17p200-h17p400
    # > status
    #   +-----------------------+-----------+-----------+-----------+-----------+
    #   | sdi                   | inactive  | submitted | success   | failure   |
    #   +-----------------------+-----------+-----------+-----------+-----------+
    #   | z22://h17p200-h17p400 | 0         | 9375      | 604       | 21        |
    #   | z22://h18p200-h18p400 | 10000     | 0         | 0         | 0         |
    #   +-----------------------+-----------+-----------+-----------+-----------+
    # > status --top
    #   +-----------------------+-----------------+-----------------+-----------+--------------------------------------------+-------------------------------------------+
    #   | sdi                   | processed poses | free poses      | top poses | poses path                                 | scores path                               |
    #   +-----------------------+-----------------+-----------------+-----------+--------------------------------------------+-------------------------------------------+
    #   | z22://h17p200-h17p400 | 9503136         | 100324          | 300000    | s3://mybucket/mylsdbase/top/a342ve.mol2.gz | s3://mybucket/mylsdbase/top/a342ve.scores |
    #   | z22://h18p200-h18p400 | 0               | 0               | 0         | N/A                                        | N/A                                       |
    #   +-----------------------+-----------------+-----------------+-----------+--------------------------------------------+-------------------------------------------+