from pydock3.files import DirBase, FileBase

from pydock3.lsd.lsd import load_lsdbase

import tarfile

# tarfile streamer- only implementing "read" because the readlines method won't be used by dock
class TarStreamer(IOBase):
    def __init__(self, tarfile : FileBase, fnfilter=None):
        tfo = tarfile.open('rb')
        self.tarfile = tarfile.Tarfile(fileobj=tfo, mode='r|gz')
        self.member = self.tarfile.extractfile(self.tarfile.next())

    def read(self, size=-1):
        b = []
        to_read = size
        while True:
            to_read = size
            b = self.member.read(to_read)
            to_read = to_read - len(b)
            if len(b) == 0 and to_read != 0:
                self.member = self.tarfile.extractfile(self.tarfile.next())

    def close(self):
        self.member.close()
        self.tarfile.close()

def run_lsd(cfg, workdir : DirBase, run : int, block : int, task_id : int):

    lsdbase = load_lsdbase(cfg)

    sdi = lsdbase.sdi(run, block)
    # gets the original sdi index (and block) associated with this task
    orig_block, orig_index = sdi.get_orig_index(task_id) 
    
    outdir = lsdbase.basedir
        .dir("lsd")
        .dir(lsdbase.block_label(orig_block))
        .dir(str(orig_index))

    outdir.create()

    entry = sdi.get(task_id)
    if entry.extension == '.tgz' or entry.extension == '.tar.gz':
        entries_fileobj = TarStreamer(entry)
    else:
        entries_fileobj = entry.open('rb')

    dockfiles = workdir.dir('dockfiles')
    DirBase.copy(lsdbase.dockfiles, dockfiles)

    indock = dockfiles.file("INDOCK")
    fix_indock(indock)

    p = sp.Popen(["dock64", indock], stdin=entries_fileobj)

    # no matter what, a file named ".lsd_status" will be created (possibly empty if the job fails)
    # this will be used as a marker for completion
    # this file is multi-purpose, it is used as a completion marker, but also contains stats about the results on success (e.g number of poses)
    FileBase.copy(workdir.file('.lsd_status'), outdir.file('.lsd_status.' + str(run)))

    r = quick_diagnose_errors(workdir)

    if r == 0:
        FileBase.copy(workdir.file("OUTDOCK"),      outdir.file("OUTDOCK." + str(run)))
        FileBase.copy(workdir.file("test.mol2.gz"), outdir.file(str(run) + ".mol2.gz"))
        if workdir.file("restart").exists():
            FileBase.copy(workdir.file("restart"),  outdir.file("restart"))
        sys.exit(0)
    sys.exit(1)

if __name__ == "__main__":

    cfgpath     = os.getenv("LSD_CFG_PATH")
    workdir     = os.getenv("LSD_WORKDIR")
    run         = os.getenv("LSD_RUN")
    block       = os.getenv("LSD_BLOCK")

    workdir = Dir(workdir)

    lsd(cfgpath, workdir, run, block)