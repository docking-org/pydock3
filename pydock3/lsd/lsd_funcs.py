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

    @staticmethod
    def cp(self, p1, p2):
        f1 = fileobj_from_path(p1)
        f2 = fileobj_from_path(p2)
        FileBase.copy(p1, p2)

    @staticmethod
    def rm(self, p):
        f = fileobj_from_path(p)
        f.rm()

class LSDBase(QBlasterBase):
    def __init__(self, cfg, filebase : ParallelJobFileBase, queue : JobQueue):
        self.filebase = filebase
        self.sdibase = filebase.dir('sdi')
        # hardcoding this here for now @ reasonable values
        self.cfg = {
            "max_queued" : 30000,
            "blocksize" : 2000
        }

    ### internal stats logic, used by jobs & status command
    # returns a string representing overall state, followed by a list of all jobs & their individual states
    # overall state can be:
    # inactive
    # running
    # 
    def get_lsd_stats(self, sdiname):
        lsd_out_dir = self.filebase.dir(sdiname)
    def get_lsd_block_stats(self, sdiname, run):
        pass
    def get_top_stats(self):
        pass

    ### queue logic
    def set_queue(self, queue_name):
        queue_obj = load_queue(queue_name)
        if not queue_obj:
            print("Could not find queue!")
            return
        else:
            queue_obj.ensure_path_is_accessible(self.filebase)
            self.queue = queue
    def get_queue(self):
        pass

    def list_queues(self):
        pass

    ### dockfiles logic
    def set_dockfiles(self, dockfiles_path, force=False):
        curr = self.get_dockfiles()
        if force or not curr:
            df_src = dirobj_from_path(dockfiles_path)
            df_dst = self.filebase.dir("dockfiles")
            DirBase.copy(df_src, df_dst)
            print("Successfully copied over dockfiles")
        else:
            rough_status = self.get_rough_lsd_status()
            if rough_status == "submitted":
                print("Error! Active jobs are currently using dockfiles. Wait for them to finish or cancel them.")
            elif rough_status in ["partial", "complete"]:
                print("Warning! You are trying to change your dockfiles but some poses have already been produced from them.")
                print("Choose one of the following actions:")
                print("1. Archive existing work, initialize a fresh run with new dockfiles")
                print("2. Destroy existing work, initialize a fresh run with new dockfiles")
                print("3. Continue with existing run, just swap out dockfiles.")
                choice = input("[choose 1, 2, or 3]: ")

                if choice == 1:
                    self.create_archive()
                    self.set_dockfiles(dockfiles_path)
                elif choice == 2:
                    self.destroy_base()
                    self.set_dockfiles(dockfiles_path)
                elif choice == 3:
                    self.set_dockfiles(dockfiles_path, force=True)

    def get_dockfiles(self):
        df = self.filebase.dir("dockfiles").file("INDOCK")
        if df.exists():
            return df.dir()
        else:
            return None

    ### archive/delete logic
    def list_archives(self):
        pass
    def create_archive(self):
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
            if lsd_status != "inactive":
                print("")
            sdi.destroy()

    def list_sdi(self):
        pass

    ### submit logic
    def __submit_sdi(self, sdi_id):
        sdi_status = get_lsd_stats(sdiname)
        sdi_status_ready = filter(sdi_status, lambda x:x[1] == "inactive")
        sdi_ready_list = [x[0] for x in sdi_status_ready]

        sdi = self.get_sdi(sdi_id)
        run = sdi.new_sdi_subset(subset=sdi_ready_list, bs=self.cfg["blocksize"])

        n_jobs = len(self.queue.list_jobs(expandarray=True))
        blk = 0
        while n_jobs < self.cfg["max_queued"]:
            blksize = len(run.block(blk))
            lsdjob = self.queue.submit("lsd", array=blksize, params={
                "LSD_BASE"  : self.filebase,
                "LSD_SDI"   : sdi_id,
                "LSD_RUN"   : run,
                "LSD_BLOCK" : blk
                "LSD_DOCKFILES" : self.get_dockfiles(),
            })
            cntjob = self.queue.submit("lsd_cnt", deps=lsdjob, params={
                "LSD_BASE"  : self.filebase
                "LSD_SDI"   : sdi_id,
                "LSD_RUN"   : run,
                "LSD_BLOCK" : blk
            })
            n_jobs += blksize + 1
            blk += 1

        self.__record_submitted(sdi_id, run)
        
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
    # qblaster lsd --base s3://mybucket/mylsdbase --creds ~/.mycreds
    # > add-sdi z22://h17p200-h17p400
    #
    # > add-sdi z22://h18p200-h18p400
    #
    # > set-dockfiles ~/mydockfiles
    #   Copying ~/mydockfiles to s3://mybucket/mylsdbase/dockfiles
    #   checksum: 123fda
    #
    # > get-dockfiles
    #   +------------------------------------+----------+
    #   | path                               | checksum |
    #   +------------------------------------+----------+
    #   | s3://mybucket/mylsdbase/dockfiles  | 123fda   |
    #   +------------------------------------+----------+
    #
    # > list-queues
    #   +------------+------------+----------+
    #   | queue name | queue type | max cap. |
    #   +------------+------------+----------+
    #   | slurm      | slurm      | 5000     |
    #   | sge        | sge        | 1000     |
    #   | aws_ue1    | AWS Batch  | 640      |
    #   +------------+------------+----------+
    #
    # > set-queue aws_ue1 # where are access keys set?
    #   Initializing aws_ue1 queue...
    #   * fetching aws config
    #   * performing one-time setup on aws account
    #   * setting up aws_ue1 environment (this may take a minute)
    #   * verifying setup
    #
    # > submit
    #   will submit 20000 jobs
    #   confirm? [y/n]: y
    #
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
    #   Creating archive @ s3://mybucket/mylsdbase/123fda.archive.zip
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
    #
    # > submit z22://h17p200-h17p400
    # > status
    #   +-----------------------+-----------+-----------+-----------+-----------+-----------+
    #   | sdi                   | inactive  | submitted | running   | success   | failure   |
    #   +-----------------------+-----------+-----------+-----------+-----------+-----------+
    #   | z22://h17p200-h17p400 | 0         | 9375      | 604       | 20        | 1         |
    #   | z22://h18p200-h18p400 | 10000     | 0         | 0         | 0         | 0         |
    #   +-----------------------+-----------+-----------+-----------+-----------+-----------+
    #
    # > status --top
    #   +-----------------------+-----------------+-----------------+-----------+--------------------------------------------+-------------------------------------------+
    #   | sdi                   | processed poses | free poses      | top poses | poses path                                 | scores path                               |
    #   +-----------------------+-----------------+-----------------+-----------+--------------------------------------------+-------------------------------------------+
    #   | z22://h17p200-h17p400 | 9503136         | 100324          | 300000    | s3://mybucket/mylsdbase/top/a342ve.mol2.gz | s3://mybucket/mylsdbase/top/a342ve.scores |
    #   | z22://h18p200-h18p400 | 0               | 0               | 0         | N/A                                        | N/A                                       |
    #   +-----------------------+-----------------+-----------------+-----------+--------------------------------------------+-------------------------------------------+