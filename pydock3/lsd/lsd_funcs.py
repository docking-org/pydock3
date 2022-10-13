import sys, os

# simple cache implementation
class SimpleFIFOCache:
    def __init__(self, maxsize=512e6):
        self.fifo = []
        self.cache = {}
        self.size = 0
        self.maxsize = maxsize
    def get(self, key):
        return self.cache.get(key)
    # allow an override length to be specified, e.g when caching a set/list/dict where len() does not accurately reflect space usage
    def set(self, key, data):
        dl = sys.getsizeof(data)
        while len(fifo) > 0 and self.size + dl > self.maxsize:
            pk = self.fifo.pop(0)
            pl = sys.getsizeof(pk)
            del self.cache[pk]
            self.size -= pl
        if len(fifo) == 0 and self.size + dl > self.maxsize:
            raise Exception('element too big for cache!')
        self.cache[key] = data
        if overridelen:
            self.lencache[key] = dl
    def clear(self):
        self.fifo = []
        self.cache = {}
        self.size = 0


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
        self.sdi_cache = SimpleFIFOCache()
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
        blk = self.sdi_cache.get(f"s{subset}b{block}")
        if not blk:
            

    def _set_cached_sdi(self, subset, block, data):
        self.sdi_cache.set(f"s{subset}b{block}", data)

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

    def subset_contains(self, subset, orig_idx):
        subset_idx_set = self.cache.get(f'set{subset}')
        if not subset_idx_set:
            subset_idx_set = set()
            for block in self.list_blocks(subset):
                for orig_idx, entry in self.list_entries(subset, block):
                    subset_idx_set.add(orig_idx)
            # we can put sets & lists in the cache too, all they need is a len() method
            self.cache.set(f'set{subset}', subset_idx_set)
        return subset_idx_set

    def basedir(self):
        return self.basedir

    def subsetdir(self, subset):
        return self.basedir().dir(f's{subset}.d')

    def blockdir(self, subset, block):
        return self.subsetdir(subset).dir(f'b{block}.sdi')

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

    def matches_signature(self, sig):
        return self.stats['sdi']['signature'].startswith(sig)

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
    STARTED=4
    SUBMITTED=5
    RUNNING=6
    UNKNOWN=7
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

class Job:
    def __init__(self, jobid, name, state):
        self.jobid = jobid
        self.name = name
        self.state = state

class SimpleSlurmQueue:
    def __init__(self, script_locations):
        self.script_locations = script_locations

    def submit(self, jobtype, jobname, args):
        script = self.script_locations[jobtype]
        subprocess.call(["sbatch", "-J", jobname, script] + args)

    def list_jobs(self):
        jobs = []
        with subprocess.Popen(["squeue", "--noheader", "-o", '%.18i %.64j, %.2t'], stdout=subprocess.PIPE) as lj:
            for line in lj.stdout.readlines():
                jobid, name, state = line.strip().split()
                if state == "R":
                    state = JobStates.RUNNING
                else:
                    state = JobStates.SUBMITTED
                jobs.append(Job(jobid, name, state))
        return jobs

    def acquire_lock(self, name, timeout):
        pass

    def release_lock(self, name):
        pass

# an extension to the SDI class which gathers and collates statistics from a compute campaign
class ComputeSDIStatsBase(SDI):

    def __init__(self, basedir, queue, outdir):
        super().__init__(basedir)
        self.cache = SimpleFIFOCache()
        self.queue = queue
        self.outdir = outdir

    # defines mapping of input idx -> output idx
    # an input idx can only be associated with one output idx (no many-many mapping)
    @abstractmethod
    def _get_output_idx(self, sdi_orig_idx):
        pass
    # does the same, but in reverse
    @abstractmethod
    def _list_input_idx(self, output_idx):
        pass
    # defines the file structure & how to fetch output statistics for a particular entry
    @abstractmethod
    def _get_output_stats_data(self, subset, output_index):
        pass

    # convenience method for listing all outputs & their inputs in a particular subset
    def get_subset_output_map(self, subset):
        out_set = self.cache.get(f'so{subset}')
        if out_set:
            return out_set

        out_set = {}
        for block in self.list_blocks(subset):
            for orig_idx, entry in self.list_entries(subset, block):
                out_idx = self._get_output_idx(orig_idx)
                if not out_set.get(out_idx):
                    out_set[out_idx] = {}
                    out_set[out_idx]['index'] = set([orig_idx])
                    out_set[out_idx]['block'] = set([block])
                else:
                    out_set[out_idx]['index'].add(orig_idx)
                    out_set[out_idx]['block'].add(block)

        self.cache.set(f'so{subset}', out_set)
        return out_set

    def _get_queue_stats(self, subset):

        output_map = self.get_subset_output_map(subset)
        output_map_rev = {}
        for o, i in output_map:
            for b in i['blocks']:
                output_map_rev[b] = o # set up reverse map
        submitted_blocks = self._get_blocks_started(subset)

        qstats = {}
        # assumed job name format:
        # {jtype}-{sdisig}-{subset}-{block}
        for job in self.queue.listjobs():
            try:
                jtype, sdisig, _subset, block = job.name.split('-')
            except: # if it doesn't match this format move on
                continue
            # filter out jobs we're not getting the stats for
            if jtype != self.jtype or not self.matches_signature(sdisig) or _subset != subset:
                continue
            else:
                qstats[output_map_rev[block]] = job.state

        for block in self._get_blocks_started(subset):
            if not qstats.get(output_map_rev[block]):
                qstats[output_map_rev[block]] = JobStates.STARTED

        return qstats

    # returns stats : dict, complete : bool
    def _get_subset_stats(self, subset):
        stats = self._load_subset_stats(subset)

        if not self._subset_is_started(subset):
            return {'all':'inactive'}
        elif self._subset_is_complete(subset):
            return stats

        sslength = self.stats['subsets'][subset]['length']
        if len(stats) == sslength:
            return stats, True

        # maps out_index -> submission state
        qstats = self._get_queue_stats(subset)
        to_save = []

        all_complete = True
        for out_index, input_set in self.get_subset_output_map(subset).items():

            if not stats.get(out_index):
                out_stats = self._get_output_stats_data(subset, out_index)
                if not out_stats:
                    all_complete = False
                    # job being in "STARTED" status here means it was submitted, but is no longer visible in the queue
                    # thus if there are no available stats for this job, it must have failed somehow
                    if   qstats.get(out_index) == JobStates.STARTED:
                        out_stats = {'state':JobStates.FAILURE}
                        to_save.append(out_index) # save this
                    # if the job is visible in the queue in any other way (SUBMITTED/RUNNING) mark it here and move on
                    elif qstats.get(out_index):
                        out_stats = {'state':qstats.get(out_index)}
                    # otherwise this job hasn't officially been submitted yet, so mark it as "STARTED"
                    else:
                        out_stats = {'state':JobStates.INACTIVE}
                else:
                    to_save.append(out_index)
                stats[out_index] = out_stats
            else:
                to_save.append(out_index)

        # only save final states, i.e success, failure
        stats_to_save = {out : stats[out] for out in to_save}
        self._save_subset_stats(subset, stats_to_save)

        if len(to_save) == sslength:
            return stats, True
        else:
            return stats, False

    def get_output_stats(self):
        allstats = {}
        for subset in self.list_subsets():
            allstats[subset] = self.get_subset_stats(subset)
        return allstats

    def _load_subset_stats(self, subset):
        statfile = self.subsetdir(subset).file('stats.json')
        return self._load_stats(statfile)

    def _save_subset_stats(self, subset, stats):
        statfile = self.subsetdir(subset).file('stats.json')
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
    def _get_subsets_started(self):
        return self._load_completion(self.basedir().file('started'))

    def _get_blocks_started(self, subset):
        return self._load_completion(self.subsetdir(subset).file('started'))

    def _mark_block_started(self, subset, block):
        self._save_completion(self.subsetdir(subset).file('started'), block)

    def _mark_subset_started(self, subset):
        self._save_completion(self.basedir().file('started'), subset)

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
        self.queue = SimpleSlurmQueue()

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
    #def set_queue(self, queue_name):
    #    queue_obj = load_queue(queue_name)
    #    if not queue_obj:
    #        print("Could not find queue!")
    #        return
    #    else:
    #        queue_obj.ensure_path_is_accessible(self.filebase)
    #        self.queue = queue
            
    #def get_queue(self):
        

    ### dockfiles logic
    def set_dockfiles(self, dockfiles_path, force=False):
        curr = self.get_dockfiles()
        if force or not curr:
            df_src = dirobj_from_path(dockfiles_path)
            df_dst = self.filebase.dir("dockfiles")
            DirBase.copy(df_src, df_dst)
            print("Successfully copied over dockfiles")
        else:

            print("Warning! You are changing dockfiles that have already been set.")
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
    ### TODO
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
        
        for sdi in self.list_sdi():
            if sdi.matches_signature(sdi_id):
                return sdi
            elif sdi.name == sdi_id:
                return sdi
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
        sdi = []
        for sdi_base in self.stats.sdi:
            sdi.append(ComputeSDIStatsBase(sdi_base))

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

        sdi = self.get_sdi(sdi_id)
        self.queue.acquire_lock(timeout=60, name=f"stats-{jtype}-{sdi_id}")

        try:
            sdi_stats = sdi.get_output_stats()

            # you are not permitted (period) to submit a job that is already submitted/running
            banned_resubmit_states = [JobStates.SUBMITTED, JobStates.RUNNING, JobStates.STARTED]
            resub_stats = {
                "inactive" : 0,
                "failed"   : 0,
                "partial"  : 0, # "partial" is an optional state for jobs like lsd that can partially complete
                "succeeded": 0
            }
            if not submitlist:
                default_resubmit_states = [JobStates.INACTIVE, JobStates.FAILED]
                sdi_ready_list = filter(lambda x: (x['state'] in default_resubmit_states) and (x['state'] not in banned_resubmit_states), sdi_stats)
                for resub_type in resub_stats:
                    resub_stats[resub_type] = len(filter(lambda x:x[1] == resub_type, sdi_ready_list))
                sdi_ready_list = [x[0] for x in sdi_ready_list]
                
            # if desired, you *may* resubmit a job that has already succeeded or is in some state other than submitted/running
            elif submitlist:
                sdi_ready_list = filter(lambda x: (flattened_status[x] not in banned_resubmit_states), submitlist)
                for resub_type in resub_stats:
                    resub_stats[resub_type] = len(filter(lambda x:flattened_status[x]==resub_type, submitlist))

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

        try:
            self.queue.acquire_lock(timeout=60, name=f"continue")

            for sdi in self.list_sdi():
                for subset in sdi.list_subsets():
                    if not sdi.subset_started(subset):
                        continue
                    for block in sdi.get_available_blocks():
                        jobname = '-'.join([self.jtype, sdi.signature, str(subset), str(block)])
                        self.queue.submit(self.jtype, jobname, sdi.blockfile(subset, block))
                        sdi._mark_block_started(subset, block)
        finally:
            self.queue.release_lock(name=f"continue")
        # psuedocode:
        # for sdi_id, subset, block in get_available_blocks_to_submit():
        #   queue.submit('lsd', array=blocksize(block), params={sdi_id, subset, block})
        #   mark_block_submitted(sdi_id, subset, block)
    
    def status(self):

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




    """ JAIL
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

    def _mark_block_complete(self, subset, block):
        self._save_completion(self.subsetdir(subset).file('complete'), block)

    def _get_blocks_complete(self, subset):
        return self._load_completion(self.subsetdir(subset).file('complete'))

    def _get_subsets_complete(self):
        return self._load_completion(self.basedir().file('complete'))

    def _mark_subset_complete(self, subset):
        self._save_completion(self.basedir().file('complete'), subset)
"""