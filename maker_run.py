########################################################################################################################
# Robert E. Jones (rej110@msstate.edu)
# GitHub Repo: https://github.com/rej110/myMAKER2-Pipeline
# usage:   maker_run.py [options]
#          run maker_run.py -h for options
# Details: This is a pipeline that automates the usage of the MAKER2 pipeline.
#          It runs RepeatModeler first followed by 4 rounds of MAKER2. This pipeline also trains ab initio gene
#          predictors and runs BUSCO between passes of MAKER.
########################################################################################################################
import os
import sys
import shutil
import argparse
import subprocess
from pathlib import Path


# TODO: Put into a config json
# System wide installed dependencies
REPEAT_MASKER_EXE = '/home/elena/miniconda3/envs/maker/bin/RepeatMasker'

# Conda installed dependencies
CONDA_MAKER_ENV = '/home/elena/miniconda3/envs/maker'

# GeneMark
GENEMARK_EXE = '/home/elena/software/gmes_linux_64_4/gmes_petap.pl'
GENEMARK_PROBUILD_EXE = '/home/elena/software/gmes_linux_64_4/probuild'

# BUSCO
RUN_BUSCO = f'{CONDA_MAKER_ENV}/bin/busco'

# database paths
BUSCO_LINEAGE = '/home/elena/data/busco_downloads/lineages/diptera_odb10'
AUGUSTUS_CONFIG = f'{CONDA_MAKER_ENV}/config'


class Paths:
    def __init__(self, pass_num):
        self.base = os.getcwd()
        
        # Gets BaseName and paths to genome, est, and protein homology data
        self.base_name, self.genome, self.est, self.prot = self.get_data_paths()

        # Sets pipeline dir
        self.pipeline_dir = f'{self.base}/pipeline'

        # RepeatModeler
        self.rep_mod_db = f'{self.pipeline_dir}/RepeatModeler'
        self.rep_mod_out = f'{self.rep_mod_db}/{self.base_name}.db-families.fa'
        # GeneMarkES
        self.gm_dir = f'{self.pipeline_dir}/GeneMarkES'
        self.gm_hmm = f'{self.gm_dir}/output/gmhmm.mod'

        # Path of current Maker run
        self.maker_dir = f'{self.pipeline_dir}/Maker_Pass{pass_num}'
        self.maker_out = f'{self.maker_dir}/{self.base_name}.maker.output'
        self.maker_log = f'{self.maker_out}/{self.base_name}_master_datastore_index.log'
        self.maker_all_gff = f'{self.maker_out}/{self.base_name}.all.gff'
        self.maker_tran = f'{self.maker_out}/{self.base_name}.all.maker.transcripts.fasta'
        self.maker_prot = f'{self.maker_out}/{self.base_name}.all.maker.proteins.fasta'

        # Path to SNAP training
        self.snap_dir = f'{self.pipeline_dir}/Snap_Training{pass_num}'
        self.snap_hmm = f'{self.snap_dir}/{self.base_name}_snap.hmm'

        # Path to files for BUSCO/AugustusTraining
        self.busco_dir = f'{self.pipeline_dir}/BUSCO_Pass{pass_num}/'

    def get_data_paths(self):
        base, genome, est, prot = '', '', '', ''
        for root, dirs, files in os.walk(self.base):
            for filename in files:
                if filename.endswith('.fna') and 'genomic' in filename:
                    base = os.path.splitext(filename)[0]  # .split('.')[0]
                    genome = os.path.join(root, filename)
                elif filename.endswith('.fna') and 'rna' in filename:
                    est = os.path.join(root, filename)
                elif filename.endswith('.faa') and 'protein' in filename:
                    prot = os.path.join(root, filename)

        return base, genome, est, prot


def execute(cmd):
    p = subprocess.run(cmd, shell=True, executable='/bin/bash')
    if p.returncode != 0:
        raise Exception('subprocess exites with non-zero exit status')


def make_dir(dir_path):
    try:
        Path(dir_path).mkdir(parents=True, exist_ok=True)
    except FileExistsError:
        shutil.rmtree(dir_path)
        Path(dir_path).mkdir(parents=True, exist_ok=True)


def run_rep_mod(threads, paths):
    make_dir(paths.rep_mod_db)
    os.chdir(paths.rep_mod_db)

    shutil.copyfile(paths.genome, f'./{paths.base_name}.db.fas')

    cmd = f'''
    BuildDatabase -name {paths.base_name}.db ./{paths.base_name}.db.fas
    RepeatModeler -pa {threads - 1} -database {paths.base_name}.db
    '''
    execute(cmd)


def get_aug_species(paths):
    aug_species = f'{AUGUSTUS_CONFIG}/species'

    make_dir(f'{aug_species}/{paths.base_name}')
    print(f'{aug_species}/{paths.base_name}')

    endings = [
        '_intron_probs.pbl',
        '_exon_probs.pbl',
        '_parameters.cfg.orig1',
        '_igenic_probs.pbl',
        '_parameters.cfg',
        '_weightmatrix.txt',
        '_metapars.cfg',
        '_metapars.utr.cfg',
        '_metapars.cgp.cfg'
    ]

    # r=root, d=directories, f = files
    for root, dirs, files in os.walk(paths.busco_dir):
        for file in files:
            for ending in endings:
                new_file = f'{aug_species}/{paths.base_name}/{paths.base_name}{ending}'

                if file.endswith(ending) is True:
                    file = os.path.join(root, file)
                    shutil.copy(file, new_file)
                    if ending is '_parameters.cfg.orig1' or ending is '_parameters.cfg':
                        old_prefix = file.split('/')[-1]

                        old_prefix = old_prefix.split('_parameters')[0]

                        f = open(new_file, 'r')
                        filedata = f.read()
                        f.close()

                        newdata = filedata.replace(old_prefix, paths.base_name)

                        f = open(new_file, 'w')
                        f.write(newdata)
                        f.close()


def run_genemark(threads):
    my_paths = Paths(pass_num=1)
    make_dir(my_paths.gm_dir)
    os.chdir(my_paths.gm_dir)

    cmd = f'''
   {GENEMARK_EXE} --ES --min_contig 10000 --cores {threads} --sequence {my_paths.genome}
    '''
    execute(cmd)


def run_snap(paths):
    # Commands to train SNAP
    # these lines go in the job script that follows if pass 1 or 2
    make_dir(paths.snap_dir)
    os.chdir(paths.snap_dir)

    cmds = [
        f'maker2zff {paths.maker_all_gff}',
        'fathom -categorize 1000 genome.ann genome.dna',
        'fathom -export 1000 -plus uni.ann uni.dna',
        'forge export.ann export.dna',
        f'hmm-assembler.pl {paths.base_name}_snap . > {paths.base_name}_snap.hmm'
    ]

    for cmd in cmds:
        execute(cmd)


def run_busco(threads, paths):
    # Final BUSCO run (likely needs to change to pass_num == 3)
    # fourth pass through maker likely not necesarry
    # Does NOT train AUGUSTUS
    make_dir(paths.busco_dir)
    os.chdir(paths.busco_dir)

    cmd = f'''
    export AUGUSTUS_CONFIG_PATH={AUGUSTUS_CONFIG}
    {RUN_BUSCO} \\
    -i {paths.maker_prot}  \\
    -l {BUSCO_LINEAGE} \\
    -o {paths.base_name}_prot \\
    -m proteins \\
    -c {threads}
    
    {RUN_BUSCO} \\
    -i {paths.maker_tran} \\
    -l {BUSCO_LINEAGE} \\
    -o {paths.base_name}_tran  \\
    -m genome \\
     -c {threads} \\
    --long
    '''
    execute(cmd)


def run_maker(new, old, threads, pass_num, alt_est):
    # Makes MAKER directory and cd's to it
    make_dir(new.maker_dir)
    os.chdir(new.maker_dir)

    # Generates and edits MAKER ctl files
    maker_ctl(pass_num=pass_num, alt_est=alt_est, old=old)

    # Runs Maker
    cmd = f'''
    maker -c {threads - 1} -fix_nucleotides maker_opts.ctl maker_bopts.ctl maker_exe.ctl
    cd {new.maker_out}
    fasta_merge -d {new.maker_log}
    gff3_merge -d {new.maker_log}
    '''

    execute(cmd)


def mod_maker_exe(old):
    exe_dict = {
        'RepeatMasker': REPEAT_MASKER_EXE,
        'probuild': GENEMARK_PROBUILD_EXE
    }
    with open('maker_exe.ctl', 'r') as infile, open('tmp.ctl', 'w') as outfile:
        for line in infile:
            line_begin = line.strip().split('=')[0]
            try:
                line_end = line.split('#')[1]
            except IndexError:
                continue

            if line_begin in exe_dict.keys():
                outfile.write(f'{line_begin}={exe_dict[line_begin]} #{line_end}')
            else:
                outfile.write(line)

    os.rename('tmp.ctl', 'maker_exe.ctl')


def maker_ctl(pass_num, alt_est, old):
    # Create the 3 MAKER2 .ctl files
    cmd = 'maker -CTL'
    execute(cmd)
    mod_maker_exe(old=old)
    
    # Initalizes data_dict
    # data_dict wil get updated as needed depening on pass_num
    data_dict = {
        'genome': old.genome,
        'rmlib': old.rep_mod_out,
        'protein': old.prot,
        'snaphmm': '',
        'gmhmm': '',
        'augustus_species': '',
        'maker_gff': '',
        'est2genome': 0,
        'protein2genome': 0,
        'keep_preds': 0
        }
    
    # Sets est evidence to est or altest depending if est evidence is from an alternate organism
    if alt_est is False:
        data_dict['est'] = old.est
    elif alt_est is True:
        data_dict['altest'] = old.est

    # Sets paths to SNAP hmm, Genemark hmm, and the AUGUSTUS species model.
    # For the "true passes" of maker. Passes 2, 3, and 4
    if pass_num in [2, 3, 4]:
        data_dict['snaphmm'] = old.snap_hmm
        data_dict['gmhmm'] = old.gm_hmm
        data_dict['augustus_species'] = old.base_name
        data_dict['maker_gff'] = old.maker_all_gff
    
    # Sets pass specific values
    if pass_num == 1:
        # sets est2genome and protein2genome depending on if that info is available
        if old.est:
            data_dict['est2genome'] = 1
        if old.prot:
            data_dict['protein2genome'] = 1
    elif pass_num == 4:
        data_dict['keep_preds'] = 1

    # Edit maker_opts.ctl - This file contatins the options for the maker run
    with open('maker_opts.ctl', 'r') as infile, open('tmp.ctl', 'w') as outfile:
        for line in infile:
            line_begin = line.strip().split('=')[0]
            try:
                line_end = line.split('#')[1]
            except IndexError:
                continue

            if line_begin in data_dict.keys():
                outfile.write(f'{line_begin}={data_dict[line_begin]} #{line_end}')
            else:
                outfile.write(line)

    # Removes original, unedited maker_opts.ctl
    os.rename('tmp.ctl','maker_opts.ctl')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Script for running MAKER2 Pipeline.',
                                     usage="maker_run.py [OPTIONS]")
    parser.add_argument('-t', '--threads', type=int, default=1,
                        help='Number of threads; default: 1')
    parser.add_argument('-p', '--passage', type=int, default=1,
                        help='Passage through MAKER2 pipeline to start from (i.e. 1, 2, 3, or 4); default: 1')
    parser.add_argument('-a', '--alt_est', default=False, action='store_true',
                        help='If est and protein data from an alternate organism; default: False')
    args = parser.parse_args()
    passage = args.passage

    for pass_num in range(passage, 5):
        new_paths = Paths(pass_num=pass_num)
        old_paths = Paths(pass_num=(pass_num - 1))

        if pass_num == 1:
            run_rep_mod(threads=args.threads, paths=new_paths)
            run_maker(threads=args.threads, pass_num=pass_num, alt_est=args.alt_est,
                      new=new_paths, old=old_paths)

            # Ab-initio training (SNAP and Augustus)
            os.chdir(old_paths.base)
            run_genemark(threads=args.threads)
            run_snap(paths=new_paths)
            run_busco(threads=args.threads, paths=new_paths)
            get_aug_species(paths=new_paths)

            os.chdir(new_paths.base)

        elif pass_num == 2 or pass_num == 3:
            run_maker(threads=args.threads, pass_num=pass_num, alt_est=args.alt_est,
                      new=new_paths, old=old_paths)

            # Ab-initio training
            run_snap(paths=new_paths)
            run_busco(threads=args.threads, paths=new_paths)
            get_aug_species(paths=new_paths)

            os.chdir(new_paths.base)

        elif pass_num == 4:
            run_maker(threads=args.threads, pass_num=pass_num, alt_est=args.alt_est,
                      new=new_paths, old=old_paths)
            run_busco(threads=args.threads, paths=new_paths)
            os.chdir(new_paths.base)
