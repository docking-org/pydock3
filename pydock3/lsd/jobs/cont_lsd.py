from pydock3.files import DirBase, FileBase

from pydock3.lsd.lsd import load_lsd

def cont_lsd(cfg, workdir, run, block, stride, priority)

    lsdbase, queue = load_lsd(cfg)
    queue.load_job_definition("lsd")
    queue.load_job_definition("cont_lsd")

    sdi = lsdbase.sdi(run, block)

    for orig_block in sdi.orig_blocks():
        blocklabel = lsdbase.block_label(block)
        blockdir   = lsdbase.dir("lsd").dir(blocklabel)
        blockdir.create()

    lsdjobname = f"lsd-r{run}-b{block}"
    cntjobname = f"cont_lsd-r{run}-b{block}"

    lsdjob = queue.submit(
        "lsd", array=sdi.blocksize, priority=priority,
        params={
            "LSD_CFG_PATH":lsdbase.cfg_path, 
            "LSD_WORKDIR":workdir,
            "LSD_RUN": run, 
            "LSD_BLOCK":block
        }
    )

    if not ((block + stride) * lsdbase.block_size > len(lsdbase.sdi)):
        queue.submit(
            "cont_lsd", deps=lsdjob, priority=priority, 
            params={
                "LSD_CFG_PATH":lsdbase.cfg_path, 
                "LSD_WORKDIR":workdir,
                "LSD_RUN": run, 
                "LSD_BLOCK":block+stride, 
                "LSD_SUBMIT_STRIDE":stride, 
                "LSD_SUBMIT_PRIO":priority
            }
        )
    return 0

if __name__ == "__main__":

    cfgpath     = os.getenv("LSD_CFG_PATH")
    workdir     = os.getenv("LSD_WORKDIR")
    run         = os.getenv("LSD_RUN")
    block       = os.getenv("LSD_BLOCK")
    stride      = os.getenv("LSD_SUBMIT_STRIDE")
    priority    = os.getenv("LSD_SUBMIT_PRIO")

    cont_lsd(cfgpath, run, block, stride, priority)