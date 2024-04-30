import argparse
import sys
import os

if __name__=="__main__":

    # get additional arguments
    parser = argparse.ArgumentParser()
    parser.add_argument('--detector', '-d', help='Detector to simulate', choices=['Par04', 'Par04_dataset2', 'ODD_dataset2', 'ODDcells_dataset2', 'FCCeeALLEGRO_dataset2', 'FCCeeCLD_dataset2', 'Par04SciPb_dataset2', 'Par04PbWO4_dataset2'], default="Par04")
    parser.add_argument('--energy', '-e',type=int, nargs="+", help='Energy to simulate (GeV)', required=True)
    parser.add_argument('--phi',type=float, help='Azimuthal angle (rad) to set if particle gun is used (= energy is non zero)',default=0)
    parser.add_argument('--theta',type=float, help='Theta angle (rad) to set if particle gun is used (= energy is non zero)',default=1.57)
    parser.add_argument('--number', '-n',type=int, default=1000, help='Number of events to simulate')
    parser.add_argument('--jobs', '-j',type=int, default=1, help='Number of jobs/runs to submit')
    parser.add_argument('--particle', '-p',type=str, default="gamma", help='Name of particle')
    parser.add_argument('--test', '-t', action='store_true', help='If chosen, only config files are created, job is not submitted')
    parser.add_argument('--flatNotPower', action='store_true', help='If chosen, and GPS is used, flat energy distribution will be generated instead of E^-1')
    parser.add_argument('--mail', '-m', help='Send mail about completed job', default="")
    args, _ = parser.parse_known_args()
    print(args)

    if len(args.energy) < 0:
    # exit if something is wrong
        sys.exit("No particle energy for simulation specified")
    if 0 in args.energy:
        if len(args.energy) > 1:
            sys.exit("Energy:=0 means use of GPS, so it cannot be combined with other energies")
        else:
            print("Use GPS instead of particle gun")
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
        if "cells" not in args.detector:
            fenergy.write('/analysis/setFileName analysis.root\n')
            fenergy.write('/param/ActivateModel DefineMeshModel\n')
            fenergy.write('/param/listModels\n')
        if energy == 0:
            fenergy.write(f'/gps/particle {args.particle}\n')
            if args.flatNotPower:
                fenergy.write('/gps/ene/type Lin\n')
                fenergy.write('/gps/ene/intercept 1\n')
            else: # use power distribution by default, E^-1
                fenergy.write('/gps/ene/type Pow\n')
                fenergy.write('/gps/ene/alpha -1\n')
            fenergy.write('/gps/ene/min 1.0 GeV\n')
            if 'FCCee' in args.detector:
                fenergy.write('/gps/ene/max 100 GeV\n')
            else:
                fenergy.write('/gps/ene/max 1.0 TeV\n')
            fenergy.write('/gps/ang/type iso\n')
            fenergy.write('/gps/ang/mintheta 0.87 rad\n')
            fenergy.write('/gps/ang/maxtheta 2.27 rad\n')
            fenergy.write('/gps/ang/minphi 0\n')
            fenergy.write('/gps/ang/maxphi 6.283 rad\n')
        fenergy.write(f'/run/beamOn {args.number}\n')
        fenergy.close()

        os.system(f'cat {input_macro}')

        # create submition script
        fsub = open('./submit.sub', 'w')
        if args.detector == 'Par04':
            if energy == 0:
                fsub.write('executable             = run_on_condor_gps.sh\n')
            else:
                fsub.write('executable             = run_on_condor_pgun.sh\n')
            detname = "Par04"
        elif args.detector == 'Par04_dataset2':
            if energy == 0:
                if args.flatNotPower:
                    fsub.write('executable             = run_on_condor_dataset2_gps_flat.sh\n')
                else:
                    fsub.write('executable             = run_on_condor_dataset2_gps.sh\n')
            else:
                fsub.write('executable             = run_on_condor_dataset2_pgun.sh\n')
            detname = "Par04"
        elif args.detector == 'ODD_dataset2':
            if energy == 0:
                if args.flatNotPower:
                    fsub.write('executable             = run_on_condor_ODD_dataset2_gps_flat.sh\n')
                else:
                    fsub.write('executable             = run_on_condor_ODD_dataset2_gps.sh\n')
            else:
                fsub.write('executable             = run_on_condor_ODD_dataset2_pgun.sh\n')
            detname = "ODD"
        elif args.detector == 'ODDcells_dataset2':
            fsub.write('executable             = run_on_condor_ODDcells_dataset2_pgun.sh\n')
            detname = "ODDcells"
        elif args.detector == 'FCCeeALLEGRO_dataset2':
            if energy == 0:
                if not args.flatNotPower:
                    fsub.write('executable             = run_on_condor_FCCeeALLEGRO_dataset2_gps.sh\n')
            else:
                fsub.write('executable             = run_on_condor_FCCeeALLEGRO_dataset2_pgun.sh\n')
            detname = "FCCeeALLEGRO"
        elif args.detector == 'FCCeeCLD_dataset2':
            if energy == 0:
                if not args.flatNotPower:
                    fsub.write('executable             = run_on_condor_FCCeeCLD_dataset2_gps.sh\n')
            else:
                fsub.write('executable             = run_on_condor_FCCeeCLD_dataset2_pgun.sh\n')
            detname = "FCCeeCLD"
        elif args.detector == 'Par04SciPb_dataset2':
            if energy == 0:
                fsub.write('executable             = run_on_condor_Par04SciPb_dataset2_gps.sh\n')
            else:
                fsub.write('executable             = run_on_condor_Par04SciPb_dataset2_pgun.sh\n')
            detname = "Par04SciPb"
        elif args.detector == 'Par04PbWO4_dataset2':
            if energy == 0:
                fsub.write('executable             = run_on_condor_Par04PbWO4_dataset2_gps.sh\n')
            else:
                fsub.write('executable             = run_on_condor_Par04PbWO4_dataset2_pgun.sh\n')
            detname = "Par04PbWO4"
        fsub.write(f'arguments             = {energy} {args.number} {args.particle} {escapedName} $(ClusterID) $(ProcId) {args.phi} {args.theta}\n')
        fsub.write(f'output                = output/job.{detname}.{escapedName}.energy{energy}.number{args.number}.id$(ClusterId).$(ProcId).out\n')
        fsub.write(f'log                   = log/job.{detname}.{escapedName}.energy{energy}.number{args.number}.id$(ClusterId).log\n')
        fsub.write(f'error                 = error/job.{detname}.{escapedName}.energy{energy}.number{args.number}.id$(ClusterId).$(ProcId).err\n')
        fsub.write('RequestCpus = 8\n')
        fsub.write('+JobFlavour = "testmatch"\n')
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
