/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package flowpro.core.solver;

import flowpro.api.Dynamics;
import flowpro.api.Equation;
import flowpro.core.Mesh;
import flowpro.core.Parameters;
import flowpro.core.Solution;
import flowpro.core.State;
import flowpro.core.parallel.AppInfo;
import flowpro.core.parallel.Domain;
import flowpro.core.parallel.Fetcher;
import flowpro.core.parallel.IpAddressReader;
import java.io.File;
import java.io.IOException;
import java.net.InetSocketAddress;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.Comparator;
import java.util.List;
import java.util.Map;
import litempi.MPIException;
import litempi.MPIMaster;
import org.apache.commons.lang3.tuple.MutablePair;

/**
 *
 * @author obublik
 */
abstract public class MasterSolver {

    public static final String APP_NAME = "app";

    abstract public Solution solve() throws MPIException, IOException;

    abstract public void saveData(Solution sol) throws IOException;

    abstract public Mesh getMesh();

    abstract public void testDynamic(double dt, int newtonIter) throws IOException;

    private static MPIMaster distributeJobs(int nSlaves, Parameters par) throws IOException, MPIException {
        IpAddressReader ipAddressReader = new IpAddressReader(par.pcFilterFile, Parameters.PC_LIST_FILE);
        List<MutablePair<String, Integer>> nodeList = ipAddressReader.getNodeList();
        Map<String, String> ip2nodeNameMap = ipAddressReader.getIp2NodeNameMap();

        Comparator<MutablePair<String, Integer>> comparator = (n1, n2) -> {
            if (n1.getValue() == n2.getValue()) {
                return 0;
            } else if (n1.getValue() > n2.getValue()) {
                return 1;
            } else {
                return -1;
            }
        };
        int maxSlavesPerNode = Collections.max(nodeList, comparator).getValue();

        String[] argArr = new String[maxSlavesPerNode];
        for (int i = 0; i < maxSlavesPerNode; i++) {
            argArr[i] = "slave " + (par.slavePort + i) + " " + par.parallelSolverType;
        }
        AppInfo appInfo = new AppInfo("FlowPro.jar", argArr, "app_");

        File zippedApp = new File("FlowPro.zip");
        Fetcher fetcher = new Fetcher(nSlaves, ip2nodeNameMap, par.fetcherPort, par.publicKeyFile);
        List<MutablePair<String, Integer>> reachedNodes = fetcher.establishConnections(nodeList);
        fetcher.fetch(zippedApp, appInfo);

        List<InetSocketAddress> inetList = new ArrayList<>();
        reachedNodes.forEach((n) -> {
            for (int i = 0; i < n.getValue(); i++) {
                inetList.add(new InetSocketAddress(n.getKey(), par.slavePort + i));
            }
        });
        InetSocketAddress[] inetAddresses = inetList.toArray(new InetSocketAddress[0]);

        return new MPIMaster(inetAddresses, par.slavePort, ip2nodeNameMap);
    }

    public enum MasterSolverType {
        ksp, schwartz;

        public static void help() {
            System.out.println("********************************");
            System.out.println("HELP for parameter spatialMethod");
            System.out.println("list of possible values:");
            System.out.println(Arrays.asList(MasterSolverType.values()));
            System.out.println("********************************");
        }
    }

    public static MasterSolver factory(String simulationPath, Mesh[] meshes, Dynamics dyn,
            Equation eqn, Parameters par, State state, Domain domain, Object lock) throws IOException, MPIException {

        MasterSolver masterSolver = null;
        if (par.parallelMode) {
            MPIMaster mpi = distributeJobs(domain.nDoms, par);
            try {
                MasterSolverType masterSolverType = MasterSolverType.valueOf(par.parallelSolverType.toLowerCase());
                switch (masterSolverType) {
                    case ksp:
                        masterSolver = new KSPSolver(mpi, simulationPath, meshes, dyn, eqn, par, state, domain, lock);
						break;
                    case schwartz:
                        masterSolver = new SchwartzImplicitSolver(mpi, simulationPath, meshes, dyn, eqn, par, state, domain, lock);
						break;
                }
            } catch (IllegalArgumentException ex) {
                MasterSolverType.help();
                throw new IOException("unknown master solver " + par.parallelSolverType.toLowerCase());
            }
            return masterSolver;
        } else {
            String solverType = ((meshes[0].getElems())[0].ti).getLocalSolverType();
            switch (solverType) {
                case "localimplicit":
                    return new LocalImplicitSolver(simulationPath, meshes, dyn, eqn, par, state, domain, lock);

                case "localexplicit":
                    return new LocalExplicitSolver(simulationPath, meshes, dyn, eqn, par, state, domain, lock);
                default:
                    throw new IOException("unknown solver " + solverType);
            }
        }
    }
}
