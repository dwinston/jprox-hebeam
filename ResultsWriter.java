package sim;

import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;

public class ResultsWriter {
	
	Simulation sim;
	String fPrefix;
	PrintWriter out;
	String filename;
	double ic;
	
	/**
	 * Real artists ship.
	 * @param sim the simulation that is done
	 * @param fPrefix  prefix of output filenames
	 */
	public ResultsWriter(Simulation sim, String fPrefix) {
		this.sim = sim;
		this.fPrefix = fPrefix;
		this.ic = (double) sim.ionCount;
	}

	public void writeResults() throws IOException {
	    String PSF = ".Ir.z"; // radial point spread function
	    String LSF = ".Ey.z"; // line spread function
	    String DSF = ".Ez"; // depth "spread" function
	    for (int k = 0; k < sim.numZs ; k++) {
	    	/* write PSF */
	    	filename = fPrefix + PSF + ((k+1)/10) + ((k+1)%10);
	    	writeSpreadFunction(k, sim.Irz, sim.DR, 1.0e24/1.60219e-19);
	    	// converted from eV/(A^3) to eV/(cm^3 C)
	    	
	        /* write LSF */
	        filename = fPrefix + LSF + ((k+1)/10) + ((k+1)%10);
	        writeSpreadFunction(k, sim.Eyz, sim.DY, 1.0e16/1.60219e-19);
	        // converted from eV/(A^2) to eV/(cm^2 C)
	    }
	    if (sim.prune) {return;}
	    /* write DSF */
        filename = fPrefix + DSF;
        writeDepthFunction(sim.Ez);
	}
	
	public void writeResults_sec() throws IOException {
		String PSF = ".Ir.sec.z"; // radial point spread function
	    String LSF = ".Ey.sec.z"; // line spread function
	    String DSF = ".Ez.sec"; // depth "spread" function    
	    /* Add primary energy contribution to
	     * Irz_sec, Eyz_sec, and Ex_sec for output */
	    for (int k = 0; k < sim.numZs;k++) {
	        for (int i = 0; i < sim.MAX_R_INDEX; i++) {
	            sim.Irz_sec[k*sim.MAX_R_INDEX + i] += sim.Irz[k*sim.MAX_R_INDEX + i];
	            sim.Eyz_sec[k*sim.MAX_R_INDEX + i] += sim.Eyz[k*sim.MAX_R_INDEX + i];
	        }
	    }
	    for (int i = 0; i < sim.MAX_Z_INDEX; i++) {
	        sim.Ez_sec[i] += sim.Ez[i];
	    }
	    
	    for (int k = 0; k < sim.numZs; k++) {
	    	/* write PSF */
	    	filename = fPrefix + PSF + ((k+1)/10) + ((k+1)%10);
	    	writeSpreadFunction(k, sim.Irz_sec, sim.DR, 1.0e24/1.60219e-19);
	    	// converted from eV/(A^3) to eV/(cm^3 C)
	    	
	        /* write LSF */
	        filename = fPrefix + LSF + ((k+1)/10) + ((k+1)%10);
	        writeSpreadFunction(k, sim.Eyz_sec, sim.DY, 1.0e16/1.60219e-19);
	        // converted from eV/(A^2) to eV/(cm^2 C)
	    }
	    if (sim.prune) {return;}
	    /* write DSF */
        filename = fPrefix + DSF;
        writeDepthFunction(sim.Ez_sec);
	}
	
	/**
	 * 
	 * @param k index for current z value in array of z values
	 * @param SF spread function being accessed
	 * @param D mesh size being used
	 * @param convFactor unit conversion factor 
	 * @throws IOException
	 */
	private void writeSpreadFunction
	(int k, double[] SF, double D, double convFactor) throws IOException {
		out = new PrintWriter(new FileWriter(filename), true);
        Simulation.out.println("Writing file: "+filename);
        int num_zero = 0; // So that more than one zero is not output
        for (int i = 0;i < sim.MAX_R_INDEX; i++) {
            if (SF[k*sim.MAX_R_INDEX + i] > 0.0) {
                if (num_zero > 1) {
                    out.format("%f %e\n",((double)i-1+0.5)*D,0.001);
                }
                num_zero = 0;
                out.format("%f %e\n",
                		((double)i+0.5)*D,
                        SF[k*sim.MAX_R_INDEX + i]/ic*convFactor);
            } else {
                if (num_zero == 0) {
                    out.format("%f %e\n",((double)i+0.5)*D,0.001);
                }
                num_zero++;
            }
        }
	}
	
	/**
	 * 
	 * @param DSF depth spread function being accessed
	 * @throws IOException
	 */
	private void writeDepthFunction(double[] DSF) throws IOException {
		out = new PrintWriter(new FileWriter(filename), true);
        Simulation.out.println("Writing file: "+filename);
	    int num_zero = 0; // So that more than one zero is not output
	    for (int i = 0; i < sim.MAX_Z_INDEX; i++) {
	    	if (sim.Ez[i] > 0.0) {
	            if (num_zero > 1)
	            	out.format("%f %e\n",((double)i-1+0.5)*sim.DZ,0.001);
	            num_zero = 0;
	            out.format("%f %e\n",
	            		((double)i+0.5)*sim.DZ,sim.Ez[i]/ic/1.60219e-19);
	            // converted from eV to eV/C
	        } else {
	            if (num_zero ==0)
	            	out.format("%f %e\n",((double)i+0.5)*sim.DZ,0.001);
	            num_zero++;
	        }
	    }
	}

}
