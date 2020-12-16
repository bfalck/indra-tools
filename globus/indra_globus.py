"""
This code is for people who know what they're doing, e.g. IDIES and SciServer admins, to
facilitate the transfer of Indra data from DataScope to FileDB or delete from FileDB using Globus.

Written by Bridget Falck, 2020.
"""

import globus_sdk
from globus_sdk import TransferClient


def get_run_num(x,y,z):
    '''Helper function to figure out raveled index from unraveled index'''
    return x*64+y*8+z

def get_xyz(run_num):
    '''Helper function to figure out unraveled index from raveled index'''
    return run_num//64, run_num//8 % 8, run_num % 8


class IndraTransfer():
    def __init__(self, token_response):
        '''
        Assume the user has logged in to globus and the jhuidies#dmztest collection, initialized their client, 
        visited the authentication url, copied the authentication code, and converted it to tokens, as
        set up in the accompanying notebook.
        e.g.:
            CLIENT_ID = 'XXXXX'
            client = globus_sdk.NativeAppAuthClient(CLIENT_ID)
            client.oauth2_start_flow()
            authorize_url = client.oauth2_get_authorize_url()
            auth_code = 'XXXXX'
            token_response = client.oauth2_exchange_code_for_tokens(auth_code)
        
        Here, set up the transfer client and jhu endpoint.
        '''
        
        globus_auth_data = token_response.by_resource_server['auth.globus.org']
        globus_transfer_data = token_response.by_resource_server['transfer.api.globus.org']

        AUTH_TOKEN = globus_auth_data['access_token']
        TRANSFER_TOKEN = globus_transfer_data['access_token']
        authorizer = globus_sdk.AccessTokenAuthorizer(TRANSFER_TOKEN)
        
        self.tc = globus_sdk.TransferClient(authorizer=authorizer)
        self.jhu_endpoint='jhuidies#dmztest'
        self.jhu_tc = self.tc.get_endpoint(self.jhu_endpoint)
        
    
    def populate_tasks(self):
        '''
        Populate lists of source (datascope) and destination (FileDB) task pairs
        
        '''
        pass
    
    def sanity_checks(self):
        '''
        Print out sanity checks, e.g., number of task pairs and a few values
        '''
        pass
    
    def submit_transfer(self, for_real = False):
        '''
        Submit tasks pairs to the transfer client
        '''
        pass
    
    def submit_delete(self, for_real = False):
        '''
        Delete files not in input skip list (e.g. keep priority snaps when deleting a full simulation)
        '''
        pass
    
    
class SimTransfer(IndraTransfer):
    '''
    Transfer one full simulation, e.g. 3_4_5, to its appropriate FileDB location. This will skip full snapshots
    that have already been transferred, and will skip sub_tab files and the FFT_DATA folder, assuming it has already
    been transferred (see the SeriesTransfer class). 
    '''
    def __init__(self, token_response, run_num):
        """
        Parameters
        ----------
        token_response : globus_sdk object
            Obtained from globus authentication code
        run_num: int
            Specifies Indra run to be transferred, from 0 to 511.
            Can be determined from X_Y_Z with indra_globus.get_run_num(X,Y,Z).
        """

        super().__init__(token_response)
        self.run_num = run_num
        self.x, self.y, self.z = get_xyz(self.run_num)
        self.fd = []
        for f in range(8,13):
            for d in range(1,4):
                self.fd.append('/sciserver/filedb{:02d}-{:02d}/cosmo/indra/'.format(f,d))
        for f in range(1,8):
            for d in range(1,4):
                self.fd.append('/sciserver/filedb{:02d}-{:02d}/cosmo/indra/'.format(f,d))

        self.dloc = '/datascope/indra{0}/{0}_{1}_{2}/'.format(self.x,self.y,self.z)
        self.floc = self.fd[self.run_num % 36]+'{}_{}_{}/'.format(self.x,self.y,self.z)
    
    def populate_tasks(self):
        '''
        Populate lists of source (datascope) and destination (FileDB) task pairs
        '''
        self.snaptasks = []
        self.subtasks = []
        self.filetasks = []

        l=[entry["name"] for entry in self.tc.operation_ls(self.jhu_tc['display_name'], path=self.floc) if entry['type'] == 'dir']
        for snapnum in range(64):
            if 'snapdir_{:03d}'.format(snapnum) in l:
                print('snapdir_{:03d} exists in {}'.format(snapnum,self.floc))
            else:
                self.snaptasks.append(('{}snapdir_{:03d}/'.format(self.dloc,snapnum),'{}snapdir_{:03d}/'.format(self.floc,snapnum)))
                for file in range(256):
                    self.subtasks.append(('{0}postproc_{1:03d}/sub_ids_{1:03d}.{2}'.format(self.dloc,snapnum,file),
                                  '{0}postproc_{1:03d}/sub_ids_{1:03d}.{2}'.format(self.floc,snapnum,file)))

        # now do text files in main directory...  powerspec_NNN.txt files, some other .txt files, and parameters-usedvalues
        l = [entry['name'] for entry in self.tc.operation_ls(self.jhu_tc['display_name'],path=self.dloc) if entry['type'] == 'file']
        for item in l:
            self.filetasks.append((self.dloc+item,self.floc+item))

    
    def sanity_checks(self):
        '''
        Print out sanity checks, e.g., number of task pairs and a few values
        '''
        print("{} snap_dir tasks, {} sub_id tasks, and {} file tasks".format(len(self.snaptasks),len(self.subtasks),len(self.filetasks)))
        print("Number of sub_id tasks should be {}".format(len(self.snaptasks)*256))
        print("First snap_dir task pair is {}".format(self.snaptasks[0]))
        print("First sub_id task pair is {}".format(self.subtasks[0]))
        print("First file task pair is {}".format(self.filetasks[0]))
        print("Last snap_dir task pair is {}".format(self.snaptasks[-1]))
        print("Last sub_id task pair is {}".format(self.subtasks[-1]))
        print("Last file task pair is {}".format(self.filetasks[-1]))
    
    
    def submit_transfer(self, for_real = False):
        '''
        Submit tasks pairs to the transfer client
        '''
        tlabel = 'Indra_{}_{}_{}'.format(self.x,self.y,self.z)
        tdata = globus_sdk.TransferData(self.tc, self.jhu_endpoint,self.jhu_endpoint,label=tlabel,sync_level="checksum")
        for task in self.snaptasks:
            tdata.add_item(task[0],task[1],recursive=True) # set recursive = True to transfer entire contents of folders
        for task in self.subtasks:
            tdata.add_item(task[0],task[1])
        for task in self.filetasks:
            tdata.add_item(task[0],task[1])

        if for_real:
            transfer_result = self.tc.submit_transfer(tdata)
            print("Transfer {} submitted!".format(tlabel))
        else:
            print("Transfer {} not submitted (set for_real=True to submit). Returning tdata.".format(tlabel))
            return tdata

    
    def submit_delete(self, skip_example = 129, for_real = False):
        '''
        Delete main dir files and snapdirs and sub_id files in non-priority snaps, which are determined
        by checking existing snapdirs in run given by skip_example (deafults to 2_0_1). This means
        that skip_example only contains priority snapshots and is not complete on FileDB.
        '''
        skipx,skipy,skipz = get_xyz(skip_example)
        skip_floc = self.fd[skip_example % 36]+'{}_{}_{}/'.format(skipx,skipy,skipz)
        
        # list snaps in skip_example: keep these
        l=[entry["name"] for entry in self.tc.operation_ls(self.jhu_tc['display_name'], path=skip_floc) if entry['type'] == 'dir']
        
        # Recursive set at upper level of DeleteData (unlike TransferData), so submit TWO Delete tasks to Globus.
        snapdata = globus_sdk.DeleteData(self.tc, self.jhu_endpoint,label='Indra_delete_{}_{}_{}_dirs'.format(self.x,self.y,self.z),recursive=True)
        filedata = globus_sdk.DeleteData(self.tc, self.jhu_endpoint,label='Indra_delete_{}_{}_{}_files'.format(self.x,self.y,self.z))
        
        for snapnum in range(64):
            if 'snapdir_{:03d}'.format(snapnum) in l:
                print('Keeping snapdir_{:03d} in {}'.format(snapnum,self.floc))
            else:
                snapdata.add_item(('{}snapdir_{:03d}/'.format(self.floc,snapnum))) # recursive = True
                for file in range(256):
                    filedata.add_item(('{0}postproc_{1:03d}/sub_ids_{1:03d}.{2}'.format(self.floc,snapnum,file)))
        l = [entry['name'] for entry in self.tc.operation_ls(self.jhu_tc['display_name'],path=self.floc) if entry['type'] == 'file']
        for item in l:
            filedata.add_item((self.floc+item))
        
        if for_real:
            print('Deleting!')
            delete_result = self.tc.submit_delete(snapdata)
            delete_result = self.tc.submit_delete(filedata)
        else:
            print('Not deleting; returning snapdata, filedata')
            print('To delete, call NAME.tc.submit_delete(snapdata) and (filedata) or set for_real=True')
            return snapdata,filedata

    

class SnapTransfer(IndraTransfer):
    '''
    Transfers full (2_0_0 to 6_7_7 (default) or 7_7_7) set of snapshots to FileDB specified by snapnum.
    Creates task pairs for sub_id files and snapdir folders and assumes sub_tab files already transferred.
    Assumes directories of specified runs already exist on FileDB (e.g. /cosmo/indra/0_0_0): if that's 
    not the case, do a SeriesTransfer first.
    '''
    def __init__(self, token_response, snapnum, include_7 = False, runfirst = None, nruns = None):
        """
        Parameters
        ----------
        token_response : globus_sdk object
            Obtained from globus authentication code
        snapnum : int
            Specify the snapshot to be transferred, from 0 to 63
        include_7: boolean (optional, defualts to False)
            If runfirst is not set, set include_7=True to transfer 2_0_0 to 7_7_7
        runfirst: int (optional, defaults to 128)
            Set this to specify (with nruns) which simulations included in transfer of given snapnum
        nruns: int (optional)
            Defaults to 320 (2_0_0 to 6_7_7); 384 (if include_7 == True); or 64 (if runfirst is set)
        """

        super().__init__(token_response)
        self.snapnum = snapnum
        self.snapstr = "%03d" % self.snapnum
        if runfirst == None:
            if include_7:
                self.runfirst = 128
                self.nruns = 384
            else: 
                self.runfirst = 128
                self.nruns = 320
        else:
            self.runfirst = runfirst
            if nruns == None:
                print('nruns not set: assuming 64')
                self.nruns = 64
            else:
                self.nruns = nruns
        print('Setting up transfer of {} runs, starting with {}, of snapnum {}'.format(self.nruns,self.runfirst,self.snapnum))

        
        self.fd = []
        for f in range(8,13):
            for d in range(1,4):
                self.fd.append('/sciserver/filedb{:02d}-{:02d}/cosmo/indra/'.format(f,d))
        for f in range(1,8):
            for d in range(1,4):
                self.fd.append('/sciserver/filedb{:02d}-{:02d}/cosmo/indra/'.format(f,d))
        self.dlocs = ['/datascope/indra{0}/{0}_{1}_{2}/'.format(i//64,i//8 % 8,i%8) for i in range(self.runfirst,self.runfirst+self.nruns)]
        self.flocs = [self.fd[i % 36]+'{}_{}_{}/'.format(i//64,i//8 % 8,i%8) for i in range(self.runfirst,self.runfirst+self.nruns)]
        
        
    def populate_tasks(self):
        '''
        Populate lists of source (datascope) and destination (FileDB) task pairs
        '''
        self.snaptasks = []
        self.subtasks = []

        for i in range(self.nruns): 
            p = self.flocs[i]
            l=[entry["name"] for entry in self.tc.operation_ls(self.jhu_tc['display_name'], path=p) if entry['type'] == 'dir']
            if 'snapdir_{}'.format(self.snapstr) in l:
                print('snapdir_{} exists in {}'.format(self.snapstr,p))
            else:
                self.snaptasks.append(('{}snapdir_{}/'.format(self.dlocs[i],self.snapstr),'{}snapdir_{}/'.format(p,self.snapstr)))
                for file in range(256):
                    self.subtasks.append(('{0}postproc_{1}/sub_ids_{1}.{2}'.format(self.dlocs[i],self.snapstr,file),
                                  '{0}postproc_{1}/sub_ids_{1}.{2}'.format(self.flocs[i],self.snapstr,file)))
    
    
    def sanity_checks(self):
        '''
        Print out sanity checks, e.g., number of task pairs and a few values
        '''
        print("{} snap_dir tasks and {} sub_id tasks".format(len(self.snaptasks),len(self.subtasks)))
        print("Number of sub_id tasks should be {}".format(len(self.snaptasks)*256))
        print("First snap_dir task pair is {}".format(self.snaptasks[0]))
        print("First sub_id task pair is {}".format(self.subtasks[0]))
        print("Last snap_dir task pair is {}".format(self.snaptasks[-1]))
        print("Last sub_id task pair is {}".format(self.subtasks[-1]))
    
    
    def submit_transfer(self, for_real = False):
        '''
        Submit tasks pairs to the transfer client
        '''
        tlabel = 'Indra_snapdir_{}'.format(self.snapstr)
        tdata = globus_sdk.TransferData(self.tc, self.jhu_endpoint,self.jhu_endpoint,
                                    label=tlabel,sync_level="checksum")
        for task in self.snaptasks:
            tdata.add_item(task[0],task[1],recursive=True) # set recursive = True to transfer entire contents of folders
        for task in self.subtasks:
            tdata.add_item(task[0],task[1])

        if for_real:
            transfer_result = self.tc.submit_transfer(tdata)
            print("Transfer {} submitted!".format(tlabel))
        else:
            print("Transfer {} not submitted (set for_real=True to submit). Returning tdata.".format(tlabel))
            return tdata


    
    def submit_delete(self, skip_list = [128,192,256,320,384], for_real = False):
        '''
        Delete snapdirs and sub_id files for one snapnum in specified range of runs, except for 
        those in skip_list (full simulations).
        '''
        skip_flocs = []
        for run_num in skip_list:
            skipx,skipy,skipz = get_xyz(run_num)
            skip_flocs.append(self.fd[run_num % 36]+'{}_{}_{}/'.format(skipx,skipy,skipz))

        # Recursive set at upper level of DeleteData (unlike TransferData), so submit TWO Delete tasks to Globus.
        snapdata = globus_sdk.DeleteData(self.tc, self.jhu_endpoint,label='Indra_delete_snap_{}_dirs'.format(self.snapstr),recursive=True)
        filedata = globus_sdk.DeleteData(self.tc, self.jhu_endpoint,label='Indra_delete_snap_{}_files'.format(self.snapstr))
        
        for i in range(self.nruns): 
            p = self.flocs[i]
            if p in skip_flocs:
                print('Keeping snapdir_{} in {}'.format(self.snapstr,p))
            else:
                snapdata.add_item(('{}snapdir_{}/'.format(p,self.snapstr)))
                for file in range(256):
                    filedata.add_item(('{0}postproc_{1}/sub_ids_{1}.{2}'.format(p,self.snapstr,file)))
        
        if for_real:
            print('Deleting!')
            delete_result = self.tc.submit_delete(snapdata)
            delete_result = self.tc.submit_delete(filedata)
        else:
            print('Not deleting; returning snapdata, filedata')
            print('To delete, call NAME.tc.submit_delete(snapdata) and (filedata) or set for_real=True')
            return snapdata,filedata

    

class SeriesTransfer(IndraTransfer):
    '''
    Transfers full series of Indra runs (e.g. all 64 7_Y_Z) to appropriate FileDB locations.
    Assumes no data exists, so first create the parent folders in /cosmo/indra/ with SeriesTransfer.make_dirs().
    Creates task pairs for FFT_DATA folders and sub_tab files, for all runs and snaps, ONLY.
    Note that for the 7 series, this should be done first, and then SnapTransfer for each desired snapnum.
    This should create 256*64(snaps)*64(sims)+64(FFT dirs) = 1048640 task pairs, < 7 million (??) maximum.
    '''
    def __init__(self, token_response, series_num):
        """
        Parameters
        ----------
        token_response : globus_sdk object
            Obtained from globus authentication code
        series_num: int
            Specifies Indra series to be transferred, from 0 to 7, e.g. 7 will transfer 7_0_0 to 7_7_7
        """

        super().__init__(token_response)
        self.series_num = series_num
        self.runfirst = self.series_num*64
        self.nruns = 64
        
        fd = []
        for f in range(8,13):
            for d in range(1,4):
                fd.append('/sciserver/filedb{:02d}-{:02d}/cosmo/indra/'.format(f,d))
        for f in range(1,8):
            for d in range(1,4):
                fd.append('/sciserver/filedb{:02d}-{:02d}/cosmo/indra/'.format(f,d))
        self.dlocs = ['/datascope/indra{0}/{0}_{1}_{2}/'.format(i//64,i//8 % 8,i%8) for i in range(self.runfirst,self.runfirst+self.nruns)]
        self.flocs = [fd[i % 36]+'{}_{}_{}/'.format(i//64,i//8 % 8,i%8) for i in range(self.runfirst,self.runfirst+self.nruns)]
        
        
    def make_dirs(self):
        '''
        Create X_Y_Z folders in appropriate FileDB locations
        '''
        for i in range(self.nruns):
            runid = self.flocs[i][-6:-1]
            p = self.flocs[i][:-6]
            l=[entry["name"] for entry in self.tc.operation_ls(self.jhu_tc['display_name'], path=p) if entry['type'] == 'dir']
            if runid in l:
                print(runid,' exists in ',p)
            else:
                self.tc.operation_mkdir(self.jhu_endpoint,self.flocs[i])
        
        
    def populate_tasks(self):
        '''
        Populate lists of source (datascope) and destination (FileDB) task pairs
        '''
        self.ffttasks = []
        self.tabtasks = []
        
        for i in range(self.nruns):
            p = self.flocs[i]
            l=[entry["name"] for entry in self.tc.operation_ls(self.jhu_tc['display_name'], path=p) if entry['type'] == 'dir']
            if 'FFT_DATA' in l:
                print('FFT_DATA exists in ',p)
            else:
                self.ffttasks.append(((self.dlocs[i]+'FFT_DATA/',self.flocs[i]+'FFT_DATA/')))
        print('Populating sub_tab task pairs (could take a while)')
        for i in range(self.nruns):
            p = self.flocs[i]
            l=[entry["name"] for entry in self.tc.operation_ls(self.jhu_tc['display_name'], path=p) if entry['type'] == 'dir']
            for snapnum in range(64):
                if 'postproc_{:03d}'.format(snapnum) not in l: # skip directory if exists, but don't print it out
                    for file in range(256):
                        self.tabtasks.append(('{0}postproc_{1:03d}/sub_tab_{1:03d}.{2}'.format(self.dlocs[i],snapnum,file),
                                  '{0}postproc_{1:03d}/sub_tab_{1:03d}.{2}'.format(self.flocs[i],snapnum,file)))
    
    
    def sanity_checks(self):
        '''
        Print out sanity checks, e.g., number of task pairs and a few values
        '''
        print("{} fft tasks and {} sub_tab tasks".format(len(self.ffttasks),len(self.tabtasks)))
        print("Number of sub_tab tasks should be {}".format(self.nruns*256*64))
        print("First fft task pair is {}".format(self.ffttasks[0]))
        print("First sub_tab task pair is {}".format(self.tabtasks[0]))
        print("Last fft task pair is {}".format(self.ffttasks[-1]))
        print("Last sub_tab task pair is {}".format(self.tabtasks[-1]))
        
    
    
    def submit_transfer(self, for_real = False):
        '''
        Submit tasks pairs to the transfer client
        '''
        tlabel = "Indra_{}".format(self.series_num)
        tdata = globus_sdk.TransferData(self.tc, self.jhu_endpoint,self.jhu_endpoint,
                            label=tlabel,sync_level="checksum")
        for task in self.ffttasks:
            tdata.add_item(task[0],task[1],recursive=True) # set recursive = True to transfer folders
        for task in self.tabtasks:
            tdata.add_item(task[0],task[1])

        if for_real:
            transfer_result = self.tc.submit_transfer(tdata)
            print("Transfer {} submitted!".format(tlabel))
        else:
            print("Transfer {} not submitted (set for_real=True to submit). Returning tdata.".format(tlabel))
            return tdata

    
    def submit_delete(self, for_real = False):
        '''
        Delete ALL contents of N_Y_Z including upper-level folders in cosmo/indra/ on FileDB.
        Assuming you have delete permissions!
        '''
        
        ddata = globus_sdk.DeleteData(self.tc, self.jhu_endpoint,label='Indra_delete_{}_series'.format(self.series_num),recursive=True)

        for i in range(self.nruns):
            ddata.add_item(self.flocs[i][:-1])
        
        if for_real:
            print('Deleting!')
            delete_result = self.tc.submit_delete(ddata)
        else:
            print('Not deleting; returning data')
            print('To delete, call NAME.tc.submit_delete(data) or set for_real=True')
            return ddata

