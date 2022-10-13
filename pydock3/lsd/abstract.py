class ComputeQueue(ABC):
    @abstractmethod
    def acquire_lock(self, lockname, timeout=None):
        pass
    @abstractmethod
    def release_lock(self, lockname):
        pass
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
    # create an archive from the current base/sdi/?subset, deactivating & removing it from the active list
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