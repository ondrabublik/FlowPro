/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package flowpro.core.solver;

import flowpro.core.Parameters;
import flowpro.core.element.Element;
import flowpro.core.element.Implicit;

/**
 *
 * @author obublik
 */
public class JacobiAssembler {

    public Element[] elems;
    private final Parameters par;
    double dt, dto;

    public JacobiAssembler(Element[] elems, Parameters par) {
        this.elems = elems;
        this.par = par;
    }

    // vytvoreni vlaken, paralelni sestaveni lokalnich matic a plneni globalni matice
    public void assemble(double dt, double dto) {  // , int newtonIter
        this.dt = dt;
        this.dto = dto;

        AssemblerThread[] assemblers = new AssemblerThread[par.nThreads];

        // vlastni vypocet, parallelni beh
        for (int v = 0; v < assemblers.length; v++) {
            assemblers[v] = new AssemblerThread(v);
            assemblers[v].start();
        }

        try {
            for (AssemblerThread assembler : assemblers) {
                assembler.join();
            }
        } catch (java.lang.InterruptedException e) {
            System.err.println(e);
            System.exit(1);
        }
    }

    private class AssemblerThread extends Thread {

        private final int id;

        AssemblerThread(int id) {
            this.id = id;
        }

        @Override
        public void run() {
            for (int i = id; i < elems.length; i += par.nThreads) {
                if (elems[i].insideComputeDomain) {
                    ((Implicit) elems[i].ti).assembleJacobiMatrix(dt, dto);
                }
            }
        }
    }
}
