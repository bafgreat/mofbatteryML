import os


def put_contents(filename, output):
    '''
    write a list object into a file
    '''
    with open(filename, 'w', encoding='utf-8') as f_obj:
        f_obj.writelines(output)
    return


def create_submit(batch, id):
    new_input = []
    new_input.append("#!/bin/bash -l\n")
    new_input.append(
        "##--------- SLURM SETTINGS TO BE ADJUSTED ---------------\n")
    new_input.append("## do not join stdout and stderr\n")
    new_input.append("#SBATCH -o job.%d.out\n" % id)
    new_input.append("#SBATCH -e job.%d.err\n" % id)
    new_input.append("\n")
    new_input.append("## name of the job\n")
    new_input.append("#SBATCH -J job-%d\n" % id)
    new_input.append("\n")
    new_input.append("## execute job from the current working directory\n")
    new_input.append("#SBATCH -D ./\n")
    new_input.append("\n")
    new_input.append("## number of nodes\n")
    new_input.append("#SBATCH -N 1\n")
    new_input.append("\n")
    new_input.append("#SBATCH --mem-per-cpu=2GB\n")
    new_input.append("\n")
    new_input.append("## time limit per node\n")
    new_input.append("#SBATCH -t 100:10:00\n")
    new_input.append("\n")
    new_input.append("## number of processors per node\n")
    new_input.append("#SBATCH -p normal\n")
    new_input.append("\n")
    new_input.append("## number of cores per processor\n")
    new_input.append("#SBATCH --ntasks-per-node=64\n")
    new_input.append("\n")
    new_input.append("#SBATCH --partition=barnard \n")
    new_input.append("\n")
    new_input.append("NTIME=500\n")
    new_input.append("\n")
    new_input.append("export OMP_NUM_THREADS=24\n")
    new_input.append("\n")
    new_input.append("source ~/.bash_profile\n")
    new_input.append("\n")
    new_input.append("module load Anaconda3/2023.07-2\n")
    new_input.append("\n")
    new_input.append("source $EBROOTANACONDA3/etc/profile.d/conda.sh\n")
    # new_input.append("\n")
    # new_input.append("source conda activate mofbatteryML\n")
    new_input.append("\n")
    new_input.append(
        "conda activate /data/horse/ws/diwo093e-mofbattryML/python=3.9\n")
    new_input.append("\n")
    # new_input.append("pip install ../../mofbatteryml/\n")
    new_input.append(
        f"python compute_complexes.py -rf ../../Battery_Result_folder -bs {batch} \n")
    base = os.getcwd()
    # submit_dir = f"submit-{id}"
    # if not os.path.exists(submit_dir):
    #     os.makedirs(submit_dir)
    # os.system(f"cp compute_complexes.py  {submit_dir}")
    # os.chdir(f"submit-{id}")
    put_contents(f"submit-{id}.sh", new_input)
    # os.chdir(base)


def create_mini_batches(data, batch_size=100):
    for i in range(0, len(data), batch_size):
        yield data[i:i+batch_size]


# Example usage:
my_list = list(range(40000))  # Example list of data

id = 0
for batch in create_mini_batches(my_list):
    create_submit(str(batch[0]) + ' ' + str(batch[-1]), id)
    id += 1
    break
    print(id)
