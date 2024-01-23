import argparse
import sys
import os

if __name__=="__main__":

    # get additional arguments
    parser = argparse.ArgumentParser()
    parser.add_argument('--detector', '-d', help='Detector to simulate', choices=['Par04'], default="Par04")
    parser.add_argument('--energy', '-e',type=int, nargs="+", help='Energy to simulate (GeV)', required=True)
    parser.add_argument('--number', '-n',type=int, default=1000, help='Number of events to simulate')
    parser.add_argument('--jobs', '-j',type=int, default=1, help='Number of jobs/runs to submit')
    parser.add_argument('--particle', '-p',type=str, default="gamma", help='Name of particle')
    parser.add_argument('--test', '-t', action='store_true', help='If chosen, only config files are created, job is not submitted')
    parser.add_argument('--mail', '-m', help='Send mail about completed job', default="")
    args, _ = parser.parse_known_args()

    if len(args.energy) < 0:
    # exit if something is wrong
        sys.exit("No particle energy for simulation specified")
    print("Number of events: ", args.number)
    known = ['e-', 'e+', 'mu-', 'mu+', 'gamma', 'geantino', 'neutron', 'proton', 'pi-', 'pi+', 'pi0']
    if not args.particle in known:
        sys.exit(f"Unknown particle <<{args.particle}>>, not found in {known}")
    escapedName = args.particle.replace('-', 'M').replace('+','P')
    print(f'Particle: {args.particle}, with the escaped name: {escapedName}')

    # create output directories if non-existent
    os.system("mkdir -p output")
    os.system("mkdir -p error")
    os.system("mkdir -p log")
    os.system("mkdir -p input")

    # submit jobs
    for energy in args.energy:
        # create input energy for simulation
        input_macro = f'./input/job.{escapedName}.{energy}GeV.N{args.number}.mac'

        fenergy = open(input_macro, 'w')
        fenergy.write('\n\n/random/setSavingFlag 1\n')
        fenergy.write('/analysis/setFileName analysis.root\n')
        fenergy.write('/param/ActivateModel DefineMeshModel\n')
        fenergy.write('/param/listModels\n')
        fenergy.write(f'/run/beamOn {args.number}\n')
        fenergy.close()

        os.system(f'cat {input_macro}')

        # create submition script
        fsub = open('./submit.sub', 'w')
        if args.detector == 'Par04':
            fsub.write('executable             = run_on_condor.sh\n')
            detname = "Par04"
        fsub.write(f'arguments             = {energy} {args.number} {args.particle} {escapedName} $(ClusterID) $(ProcId)\n')
        fsub.write(f'output                = output/job.{detname}.{escapedName}.energy{energy}.number{args.number}.id$(ClusterId).$(ProcId).out\n')
        fsub.write(f'log                   = log/job.{detname}.{escapedName}.energy{energy}.number{args.number}.id$(ClusterId).log\n')
        fsub.write(f'error                 = error/job.{detname}.{escapedName}.energy{energy}.number{args.number}.id$(ClusterId).$(ProcId).err\n')
        fsub.write('RequestCpus = 8\n')
        fsub.write('+JobFlavour = "nextweek"\n')
        fsub.write(f'transfer_input_files    = {input_macro}\n')
        fsub.write(f'transfer_output_files   = ""\n')

        if args.mail != "":
            fsub.write(f'notify_user = {args.mail}\nnotification = always\n')
        fsub.write(f'queue {args.jobs}\n')
        fsub.close()

        # check what is submitted
        os.system("cat submit.sub")

        # submit job to condor
        if not args.test:
            os.system("condor_submit submit.sub")
