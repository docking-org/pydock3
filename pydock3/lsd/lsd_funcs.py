# standalone SDI implementation w/subsets, blocks, metadata, caching, etc.
# clarification of terms:
#
# sdi:      standard database index- "a list of files" or even just simply "a list", the fancy name connotes the context this list is used in
#           an "sdi" is a list that is given as the input for a distributed program, most often in a large compute cluster environment.
#
# subset:   a subset of the sdi list- subset #0 will always be the sdi list itself
#           it is useful to create subsets in the sdi context, as this allows us to create distinct groupings of list items that can be treated separately
#           for example, if in a large distributed job of 10000, 9943 succeed and 57 fail, we would want to create a list of just the 57 items that failed for re-submission (whilst still being treated as a part of the parent sdi)
#
# blocks:   fixed-length shards of a subset. these form the "atoms" of the SDI construct
#           the reasoning behind blocks is practical- if job 999995 of 999999 wants to find the list item associated with it's job index, it needs to read the entry from a source file
#           if this source file contains all 999999 lines, this job will need to read almost every single one of them, which is computationally (and monetarily) expensive
#           thus, we break the source file into predictable fixed-size blocks
#           this also aligns with job submission patterns, where jobs are submitted in fixed-size arrays with a maximum number of tasks per job array
#
# metadata: simply overall statistics about the sdi construct. Number of subsets, blocks, how many entries per a block, how many blocks per a subset, etc. 
#           stored in a single file; $SDI_BASE/sdistats
#
class SDI:
    def __init__(self, basedir):
        self.basedir = basedir
        self.blocksize = 2000
        self.sdi_cache = {}
        self.sdi_cache_fifo = []
        self.cache_limit = 1000 # hold up to 1000 blocks in memory

    def create(self, sourcefile, name=None):
        self.basedir().create()
        orig_file = self.basedir.file('orig')
        if orig_file.exists():
            raise Exception(f'Unable to create sdi - it seems that it has already been created before!\nRemove the SDI before creating again.')

        source_fo = sourcefile.open('r')
        orig_fo   = orig_file.open('w')
        self.stats = {
            'sdi' : {
                'length' : 0,
                'signature' : None,
                'name' : None
            }
            'subsets' : []
        }
        with source_fo, orig_fo:
            source_data = source_fo.read()
            source_lines = source_data.split('\n')
            orig_fo.write(source_data)

            self.stats['sdi']['length'] = len(source_lines)
            self.stats['sdi']['signature'] = hashlib.sha256(source_data.encode('utf-8')).hexdigest()
            self.stats['sdi']['name'] = name or sourcefile.basename()
        
        self.save_stats()
        # add default subset 0
        self.add_subset(list(range(len(source_lines))))

    def _get_cached_sdi(self, subset, block):
        return self.sdi_cache.get(f"s{subset}b{block}")

    def _set_cached_sdi(self, subset, block, data):
        # pop oldest entry when we reach the limit
        if len(self.sdi_cache) + 1 > self.cache_limit:
            pop_k = self.sdi_cache_fifo.pop()
            self.sdi_cache.pop(pop_k)
        self.sdi_cache[f"s{subset}b{block}"] = data
        self.sdi_cache_fifo.append()

    def get_entry(self, subset, block, index) -> int, str, DirBase:
        sdi_lines = self._get_cached_sdi(subset, block)
        if not sdi_lines:
            blockfile = self.sdidir_base.dir('s{subset}.d').file('b{block}.sdi')
            if not blockfile.exists():
                raise Exception("tried to access non-existent sdi block!")
            with blockfile.open('r') as bf:
                sdi_lines = bf.readlines()
            self.set_cached_sdi(subset, block, sdi_lines)
        orig_idx, entry_val = sdi_lines[index].split('\t')
        return orig_idx, entry_val

    def list_entries(self, subset, block) -> list[str]:
        sdi_lines = self._get_cached_sdi(subset, block)
        return sdi_lines

    def list_blocks(self, subset):
        ss_list = self.list_subsets()
        return [i for i, e in enumerate(self.stats['subsets'][subset]['blocks'])]

    def list_subsets(self):
        return sorted([r['ssid'] for r in self.stats['subsets']])
            
    def add_subset(self, subset_indices : list[int]):
        subset_id = len(self.stats['subsets'])

        new_record = {}
        new_record['ssid'] = subset_id
        new_record['length'] = len(subset_indices)
        new_record['blocks'] = []

        orig_file = self.basedir().file('orig')
        new_lines = []
        all_lines = []

        with orig_file.open('r') as of:
            all_lines = of.readlines()

        for idx in subset_indices:
            new_lines.append(f"{idx}\t{all_lines[idx]}")

        block = 0
        i = 0
        currfile = self.blockdir(subset, block).file('sdi')
        currhndl = currfile.open('w')
        while i < len(new_lines):
            currhndl.write(new_lines[i])
            i += 1
            if i % self.blocksize == 0:
                new_record['blocks'].append(blocksize)
                block += 1
                currhndl.close()
                currfile = self.blockdir(subset, block).file('sdi')
                currhndl = currfile.open('w')

        currhndl.close()
        new_record['blocks'].append(i % self.blocksize)
        new_record['nblocks'] = len(new_record['blocks'])

        self.stats['subsets'].append(new_record)
        self.save_stats()
        return subset

    def basedir(self):
        return self.basedir

    def subsetdir(self, subset):
        return self.basedir().dir('s{subset}.d')

    def blockdir(self, subset, block):
        return self.subsetdir(subset).dir('b{block}.sdi')

    @property
    def stats(self):
        if not self._stats:
            self._stats = self._get_records()
            return self._stats
        else:
            return self._stats

    @stats.setter
    def stats(self, st):
        self._stats = st

    def save_stats(self):
        self._save_records(self._stats)

    # sdi records mockup:
    """
    {
        'sdi' : {
            'length' : int(),
            'signature' : str(),
            'name' : str()
        },

        'subsets' : [
            {
                'ssid' : int(),
                'length' : int(),
                'blocksize' : int(),
                # these could technically be calculated from the above terms, but who cares
                'nblocks' : int(),
                'lastblocksize' : int()
        ]
    }
    """
    def _get_records(self):
        records_file = self.basedir().file('sdistats')
        return self._load_stats(records_file)

    def _save_records(self, stats):
        records_file = self.basedir().file('sdistats')
        self._save_stats(records_file, stats)

    # cloning this code over here so SDI can keep records and such...
    def _load_stats(self, statfile):
        if not statfile.exists():
            return {}
        with statfile.open('r') as stat_fo:
            stats = json.loads(stat_fo.read())
            return stats

    def _save_stats(self, statfile, stats):
        with statfile.open('w') as stats_fo:
            json.dump(stats, stats_fo)

class JobStates(Enum):
    
    SUCCESS=0
    FAILURE=1
    PARTIAL=2
    INACTIVE=3
    SUBMITTED=4
    RUNNING=5
    UNKNOWN=6
    def from_string(s):
        s_dict = {
            "INACTIVE":JobStates.INACTIVE,
            "SUBMITTED":JobStates.SUBMITTED,
            "RUNNING":JobStates.RUNNING,
            "SUCCESS":JobStates.SUCCESS,
            "FAILURE":JobStates.FAILURE,
            "PARTIAL":JobStates.PARTIAL
        }
        if not s_dict.get(s):
            return JobStates.UNKNOWN
        else:
            return s_dict[s]

# an extension to the SDI class which gathers and collates statistics from a compute campaign
# different implementations 
class ComputeSDIStatsBase:

    def __init__(self, sdi, queue, outdir):
        pass

    def _get_output_stats(self, output_index):


    def _get_queue_stats(self, subsetfilter=None):
        sdim = {}
        qstats = {}

        for job in self.queue.listjobs():
            jtype, sdi_id, subset, block = j.name.split('-')
            subset = int(subset)
            block = int(block)
            if jtype != self.jtype:
                continue
            if sdifilter and sdi_id != sdifilter:
                continue
            if subsetfilter and subset != subsetfilter:
                continue
            sdi = sdim.get(sdi_id)
            if not sdi:
                sdi = self.get_sdi(sdi_id)
                sdim[sdi_id] = sdi
                qstats[sdi_id] = {}
            orig_idx = sdi.get_orig_idx(subset, block, j.task_id)
            qstats[sdi_id][orig_idx] = j.state 

        if sdifilter:
            return qstats.get(sdifilter) or {}
        return qstats

    def _get_all_stats(self):
        qstats = self._get_queue_stats()
        stats_all = []
        for sdi in self.list_sdi():
            qstats_sdi = qstats.get(sdi.id) or {}
            stats_all.append(self._get_sdi_stats(sdi, qstats=qstats_sdi))
        return stats_all

    def _get_sdi_stats(self, sdi, qstats=None):
        stats = {}
        if not qstats:
            qstats = self._get_queue_stats(sdifilter=sdi.id)

        started_subsets = self._get_subsets_started(sdi)
        complete_subsets = self._get_subsets_complete(sdi)

        for subset in sdi.list_subsets():
            complete = subset in complete_subsets
            started  = subset in started_subsets
            if started:
                stats.update(self._get_subset_stats(sdi, subset, qstats=qstats, complete=complete))
        return stats

    def _get_subset_stats(self, sdi, subset, qstats=None, complete=False):
        if not qstats:
            qstats = self._get_queue_stats(sdi, sdifilter=sdi.id, subsetfilter=subset)

        stats = self._load_subset_stats(sdi, subset)
        if complete:
            return stats

        ss_complete = True
        started_blocks = self._get_blocks_started(sdi, subset)
        complete_blocks = self._get_blocks_complete(sdi, subset)

        for block in sdi.list_blocks(subset):
            if (block in complete_blocks):
                continue
            if not (block in started_blocks):
                for orig_idx, entry in sdi.list_entries(subset, block):
                    qstats[orig_idx] = JobStates.SUBMITTED
            bk_complete = True
            for orig_idx, entry in sdi.list_entries(subset, block):
                if qstats.get(orig_idx):
                    bk_complete = False
                    continue
                elif stats.get(orig_idx):
                    continue
                else:
                    stats.update(self._get_entry_stats(sdi, subset, orig_idx))
            if not bk_complete:
                ss_complete = False
            else:
                self._mark_block_complete(sdi, subset, block)

        if ss_complete:
            self._mark_subset_complete(sdi, subset)
        self._save_subset_stats(sdi, subset, stats)

        return stats

    def _get_entry_stats(self, sdi, subset, orig_index):
        outdir = sdi.get_entry_outdir(index)
        statsf = outdir.file(f'{subset}-stats.json')
        if not statsf.exists():
            return JobStates.INACTIVE, {}
        with statsf.open('r') as statsf_obj:
            stats = json.loads(statsf_obj.read())
        return stats

    def _load_subset_stats(self, sdi, subset):
        statfile = sdi.subsetdir(subset).file('stats.json')
        return self._load_stats(statfile)

    def _save_subset_stats(self, sdi, subset, stats):
        statfile = sdi.subsetdir(subset).file('stats.json')
        return self._save_stats(statfile, stats)

    def _load_stats(self, statfile):
        if not statfile.exists():
            return {}
        with statfile.open('r') as stat_fo:
            stats = json.loads(stat_fo.read())
            return stats

    def _save_stats(self, statfile, stats):
        with statfile.open('w') as stats_fo:
            json.dump(stats, stats_fo)

    # functions for setting various completion markers
    # all of these assume a lock has been acquired for exclusive access to relevant files
    def _get_subsets_started(self, sdi):
        return self._load_completion(sdi.basedir().file('started'))

    def _get_subsets_complete(self, sdi):
        return self._load_completion(sdi.basedir().file('complete'))

    def _get_blocks_started(self, sdi, subset):
        return self._load_completion(sdi.subsetdir(subset).file('started'))

    def _get_blocks_complete(self, sdi, subset):
        return self._load_completion(sdi.subsetdir(subset).file('complete'))

    def _mark_block_started(self, sdi, subset, block):
        self._save_completion(sdi.subsetdir(subset).file('started'), block)

    def _mark_block_complete(self, sdi, subset, block):
        self._save_completion(sdi.subsetdir(subset).file('complete'), block)

    def _mark_subset_started(self, sdi, subset):
        self._save_completion(sdi.basedir().file('started'), subset)

    def _mark_subset_complete(self, sdi, subset):
        self._save_completion(sdi.basedir().file('complete'), subset)

    def _load_completion(self, cfile):
        with cfile.open('r') as cfile_o:
            data = cfile_o.read().split('\n')
            return [int(d) for d in data]

    def _save_completion(self, cfile, citem):
        with cfile.open('a') as cfile_o:
            cfile_o.write(f'{citem}\n')

class LSDBase(QBlasterBase):
    def __init__(self, cfg, filebase : ParallelJobFileBase, queue : JobQueue):
        self.filebase = filebase
        self.sdibase = filebase.dir('sdi')
        self.jtype = "lsd"

    

    ### internal stats logic, used by jobs & status command
    # returns a string representing overall state, followed by a list of all jobs & their individual states
    # overall state can be:
    # state name  | allowed to submit?
    # inactive    | yes
    # running     | no
    # submitted   | no
    # failed      | yes
    # succeeded   | optional

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
                choice = input("[choose from 1, 2, or 3]: ")

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

        sdi_record = self.sdibase.file('sdilist')
        sdi_stats = self._load_stats(sdi_record)

        new_sdi = SDI.from_path(sdi_file)
        if self.get_sdi(new_sdi.signature):
            print("One of two things just happened-")
            print("Either you just tried to upload an SDI file that already exists, OR")
            print("You just experienced a thermodynamic miracle in the form of a sha256 hash collision, in which case you should drop what you're doing and go buy some lottery tickets")
            print(f"hash={new_sdi.signature}")
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
    ### funny enough, doesn't *actually* submit jobs to the queue, it just marks them as such
    ### a few words of explanation:
    ###
    ### often there is a limit on how many jobs are allowed to exist in the pyhsical queue at the same time
    ### for example, SLURM has a configuration-defined upper limit on the number of jobs that can be queued up at any given time
    ### this limit may be surpassed many times over by a big enough job, thus we cannot fit everything in the queue at once
    ### there are two approaches I know of to circumvent this:
    ###
    ### in approach 1 you maintain a constantly running "daemon" process that monitors the queue & submits more jobs when possible
    ### this approach is simple and gets results, but is not without its pitfalls
    ### the daemon process itself is vulnerable to interruption, and requires a dedicated CPU to run the polling and submission
    ###
    ### in approach 2 (the one we take here) jobs are instead marked as submitted in a log
    ### the actual submission occurs in a "continuation" function
    ### this function submits as many blocks as possible from jobs marked as submitted until the queue is full
    ### blocks that are actually submitted get marked as such- this is important for determining the trajectory of a job in transit
    ### we can submit continuation jobs as dependencies to compute blocks such that the queue is refreshed often without requiring a dedicated process
    ### this also means we need to have a lock accessible from the queue
    def __submit_sdi(self, sdi_id, dockfilesname, submitlist=None):

        self.queue.acquire_lock(timeout=60, name=f"stats-{jtype}-{sdi_id}")

        try:
            sdi_blocks_status, sdi_flattened_status = self.get_lsd_stats(sdi_id)

            # you are not permitted (period) to submit a job that is already submitted/running
            banned_resubmit_states = ["submitted", "running"]
            resub_stats = {
                "inactive" : 0,
                "failed"   : 0,
                "partial"  : 0, # "partial" is an optional state for jobs like lsd that can partially complete
                "succeeded": 0
            }
            if not submitlist:
                default_resubmit_states = ["inactive", "failed"]
                sdi_ready_list = filter(lambda x: (x[1] in default_resubmit_states) and (x[1] not in banned_resubmit_states), sdi_flattened_status)
                for resub_type in resub_stats:
                    resub_stats[resub_type] = len(filter(lambda x:x[1] == resub_type, sdi_ready_list))
                sdi_ready_list = [x[0] for x in sdi_ready_list]
                
            # if desired, you *may* resubmit a job that has already succeeded or is in some state other than submitted/running
            elif submitlist:
                sdi_ready_list = filter(lambda x: (flattened_status[x] not in banned_resubmit_states), submitlist)
                for resub_type in resub_stats:
                    resub_stats[resub_type] = len(filter(lambda x:flattened_status[x]==resub_type, submitlist))

            sdi = self.get_sdi(sdi_id)
            subset = sdi.add_subset(subset_indices=sdi_ready_list)

            # no actual submission happens yet- we merely mark it as such
            # this way our continuation submission has a record of all jobs to submit blocks for
            self._mark_subset_started(sdi, subset)

        finally:

            self.queue.release_lock(name=f"submit-{jtype}-{sdi_id}")

        self.continue_submission()

    def submit(self, sdiname=None, confirm=False):

        if sdiname is None:
            for sdiname in self.listsdi():
                self.__submit_sdi(sdiname)
        else:
            self.__submit_sdi(sdiname)

        self.continue_submission()

    # continuation function- very important, does the actual submission to queue. "submit" just gathers statistics & creates new subsets/submit markers for continuation to use
    # continue does not operate on a single sdi/subset, it operates on the submission queue as a whole
    # ergo all sdi/subsets that have been submitted are *actually* submitted through this function
    # this allows all jobs to share the same queue, while not requiring a dedicated process that constantly exists and monitors for free space in the queue per sdi
    def continue_submission(self):

        self.queue.acquire_lock(timeout=60, name=f"continue")

        # psuedocode:
        # for sdi_id, subset, block in get_available_blocks_to_submit():
        #   queue.submit('lsd', array=blocksize(block), params={sdi_id, subset, block})
        #   mark_block_submitted(sdi_id, subset, block)
    
    def status(self, sdiname=None):

        pass

    ###################
    # example session #
    ###################

    # qblaster ligands --base /nfs/home/mystuff

    # qblaster ifp

    # qblaster ecfp

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
    # > submit z22://h17p200-h17p400
    # > status
    #   +-----------------------+-----------+-----------+-----------+-----------+-----------+
    #   | sdi                   | inactive  | submitted | running   | success   | failure   | poses     | total time  | poses data size | total data processed |
    #   +-----------------------+-----------+-----------+-----------+-----------+-----------+
    #   | z22://h17p200-h17p400 | 0         | 9375      | 604       | 20        | 1         |
    #   | z22://h18p200-h18p400 | 10000     | 0         | 0         | 0         | 0         |
    #   +-----------------------+-----------+-----------+-----------+-----------+-----------+
    # > status --top
    #   +-----------------------+-----------------+-----------------+-----------+--------------------------------------------+-------------------------------------------+
    #   | sdi                   | harvested poses | free poses      | top poses | poses path                                 | scores path                               |
    #   +-----------------------+-----------------+-----------------+-----------+--------------------------------------------+-------------------------------------------+
    #   | z22://h17p200-h17p400 | 9503136         | 100324          | 300000    | s3://mybucket/mylsdbase/top/a342ve.mol2.gz | s3://mybucket/mylsdbase/top/a342ve.scores |
    #   | z22://h18p200-h18p400 | 0               | 0               | 0         | N/A                                        | N/A                                       |
    #   +-----------------------+-----------------+-----------------+-----------+--------------------------------------------+-------------------------------------------+
    # > cp 
    # > rm
    # > 