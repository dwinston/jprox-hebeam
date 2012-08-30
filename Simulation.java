package sim;

import java.io.BufferedReader;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayDeque;
import java.util.Deque;
import java.util.Random;
import java.util.Scanner;
import static java.lang.Math.abs;
import static java.lang.Math.max;
import static java.lang.Math.sqrt;
import static java.lang.Math.log;
import static java.lang.Math.pow;
import static sim.Standard.PI;
import static sim.Standard.TWOPI;
import static sim.Standard.PIOVERTWO;
import static sim.Standard.EPS;
import static sim.Standard.LAMBDA_D;

/**
 * This is the main class for running a simulation. Simulation parameters are specified within.
 * @author Donny
 *
 */
public class Simulation {
	Scanner input;
	static PrintWriter out = new PrintWriter(System.out, true);
	final int MAX_R_INDEX = 16000;
	final int MAX_Z_INDEX = 16000;
	final double E_MIN = 50;
	final int[] z = new int[MAX_Z_INDEX+1];
	final double[] Ez = new double[MAX_Z_INDEX+1];
	final double[] Ez_sec = new double[MAX_Z_INDEX+1];
	static Random rand;
	
	// Prune paths?
	boolean prune;

	Sample sample;
	
	double DR, DY, DZ, zMax;
	int numZs;
	double[] Irz, Eyz, Irz_sec, Eyz_sec;
	int ionCount;
	double initialStoppingPower;
/*	double beamEnergy, beamDiameter;
	int bsiCount = 0; // count backscattered ions
*/	Deque<Particle> queue;
	
	/**
	 * @param args
	 */
	/**
	 * @param args
	 */
	public static void main(String[] args) {
		if (args.length != 2) {
			System.err.println("Enter exactly two filenames as arguments: "
					+" e.g. 'sim_HSQ45.Si' and 'sim_COLLISION-fancy.txt'.");
			return;
		}
		out.println(args[0]+" "+args[1]);
		Simulation sim = new Simulation();
		try {
			sim.input = new Scanner(
					new BufferedReader(new FileReader(args[0])));
		} catch (FileNotFoundException e) {
			System.err.println("Could not open input file.");
			return;
		}
		rand = new Random();
		sim.configure(sim.input); out.println("");
		sim.getSampleInfo(sim.input); out.println("");
		sim.getMeshInfo(sim.input); out.println("");
		sim.getBeamInfo(sim.input); out.println("");
		sim.input.close();
		
		try {
			sim.input = new Scanner(
					new BufferedReader(new FileReader(args[1])));
		} catch (FileNotFoundException e) {
			System.err.println("Could not open input file.");
			return;
		}
		sim.queue = new ArrayDeque<Particle>();
		int pNum = 0;
		// Need to grab first mention of first ion number
		// because processOneIon() grabs first mention of
		// subsequent ion numbers.
		sim.input.nextInt();
		while (sim.input.hasNext()) {
			if (pNum % 1000 == 0) {
				out.println(pNum+" ions done");
			}
			pNum += 1;
			// possibly add SEs to sim.queue
			sim.processOneIon(pNum);
			Particle joe = sim.getNextJoe();
			// deplete sim.queue
			while (joe != null) {
				sim.calcTrajectory(joe);
				joe = sim.getNextJoe();
			}
		}
		assert (sim.ionCount == pNum);
		sim.input.close();
		
		ResultsWriter writer = new ResultsWriter(sim, args[0]);
		Simulation.out.println("Root output filename: "+args[0]);
		try {
			writer.writeResults();
			if (Particle.maxRank > 0) {writer.writeResults_sec();}
		} catch (IOException e) {
			out.println("Unable to access output file.");
		}
/*		out.println("BSI: eta="
				+(((double) sim.bsiCount)/sim.ionCount));*/
	}

	private void configure(Scanner input) {
		prune = (input.nextInt() != 0) ? true: false;
		out.println("Prune paths? (1->true, 0->false): "+prune);
		Particle.maxRank = input.nextInt();
		out.println("How many levels of secondaries? "
				+ "(0 means no secondaries): " + Particle.maxRank);
		if (Particle.maxRank > 0) {
			Particle.cutoffE = input.nextDouble();
			out.println("Cutoff energy (eV): "+Particle.cutoffE);
		}
	}
	
	private void getSampleInfo(Scanner input) {
		int numLayers = input.nextInt();
		out.println("Number of layers: "+numLayers);
		sample = new Sample(numLayers);
		for (int i = 0; i < numLayers; i++) {
			out.println("Layer "+(i+1)+":");
			Layer layer = new Layer();
			layer.thickness = input.nextDouble();
			out.println("  Thickness (A): "+layer.thickness);
			layer.numSpecies = input.nextInt();
			out.println("  Number of atom species in layer: "
					+ layer.numSpecies);
			layer.species = new Species[layer.numSpecies];
			for (int j = 0; j < layer.numSpecies; j++) {
				out.println("  Species "+(j+1)+":");
				Species species = new Species();
				species.atomicNumber = input.nextDouble();
				out.println("    Atomic number: "+species.atomicNumber);
				species.atomicWeight = input.nextDouble();
				out.println("    Atomic weight: "+species.atomicWeight);
				species.speciesDensity = input.nextDouble();
				out.println("    Density of species in layer (g/cm^3): "
						+ species.speciesDensity);
				species.update();
				out.println("    Species excitation energy (eV): "
						+ species.excitationEnergy);
				layer.species[j] = species;
			}
			layer.update();
			out.format("  Mean excitation energy (eV): %.1f%n", 
					layer.meanExcitationEnergy);
			sample.layer[i] = layer;
		}
		sample.update();
	}

	private void getMeshInfo(Scanner input) {
		DR = input.nextDouble();
		DY = DR;
		out.println("dr, dy mesh size (A): "+DR);
		DZ = input.nextDouble();
		out.println("dz mesh size (A): "+DZ);
		numZs = input.nextInt();
		out.println("Number of z values for output: "+numZs);
		Irz = new double[numZs * MAX_R_INDEX]; // could be + numZs. Eh.
		Eyz = new double[numZs * MAX_R_INDEX];
		Irz_sec = new double[numZs * MAX_R_INDEX];
		Eyz_sec = new double[numZs * MAX_R_INDEX];
		for (int i=0; i < (numZs * MAX_R_INDEX); i++) {
			Irz[i] = 0; Eyz[i] = 0; Irz_sec[i] = 0; Eyz_sec[i] = 0;
		}
		for (int i=0; i <= MAX_Z_INDEX; i++) {
			z[i] = 0; Ez[i] = 0; Ez_sec[i] = 0;
		}
		zMax = 0;
		for (int i = 0; i < numZs; i++) {
			double zValue = input.nextDouble();
			out.format("  z[%d]: %.1f%n", i, zValue);
			int zIndex = (int) (zValue/DZ + 0.5);
			if (zIndex < MAX_Z_INDEX) {
				// +1 so can check as a boolean
				z[zIndex] = i+1;
				zMax = Math.max(zMax, zValue);
			} else {
				System.err.println(""+zValue+" out of mesh. "
						+ "Please input another value.");
			}
		}
	}
	
	private void getBeamInfo(Scanner input) {
		//beamEnergy = input.nextDouble();
		//out.println("Beam energy (eV): "+beamEnergy);
		//beamDiameter = input.nextDouble();
		//out.println("Beam diameter (A): "+beamDiameter);
		ionCount = input.nextInt();
		out.println("Number of ions: "+ionCount);
		initialStoppingPower = input.nextDouble();
		out.println("Initial elec. stopping power in top layer (eV/A): "+initialStoppingPower);
	}
	
/*	private Particle getRandomJoe() {
		Particle joe = new Electron();
		joe.energy = beamEnergy;
		double phi = 2 * Standard.PI * rand.nextDouble();
		// Gaussian random variable with 
		// standard deviation = (beamDiameter/2)
		double R = (beamDiameter/2) * Math.abs(rand.nextGaussian());
		joe.rNew = new Coord(R * Math.sin(phi),
				R * Math.cos(phi),
				0.001); // offset to avoid division by zero
		// propagate normal to surface to first scattering event
		joe.phi = joe.theta = 0;
		joe.r = new Coord(joe.rNew.x, joe.rNew.y, 0);
		return joe;
	}*/
	
	private void processOneIon(int ionNum) {
		Coord r = new Coord(0,0,0); // initial ion position
		double stoppingPower_r = initialStoppingPower;
		input.nextDouble(); // ion energy (in keV), unused for now
		double z = input.nextDouble();
		double y = input.nextDouble();
		double x = input.nextDouble();
		Coord rNew = new Coord(x,y,z);
		double stoppingPower_rNew = input.nextDouble(); // in eV/Angstrom
		double avgSP = (stoppingPower_r + stoppingPower_rNew) / 2;
		if (Particle.maxRank > 0) {
			generateSEs(r,rNew,avgSP);
		} else {
			dissipateSwissly(r,rNew,avgSP);
		}
		r = rNew;
		stoppingPower_r = stoppingPower_rNew;
		while ((input.hasNext()) && (input.nextInt() == ionNum)) {
			input.nextDouble(); // ion energy (in keV), unused for now
			z = input.nextDouble();
			y = input.nextDouble();
			x = input.nextDouble();
			rNew = new Coord(x,y,z);
			stoppingPower_rNew = input.nextDouble();
			avgSP = (stoppingPower_r + stoppingPower_rNew) / 2;
			if (Particle.maxRank > 0) {
				generateSEs(r,rNew,avgSP);
			} else {
				dissipateSwissly(r,rNew,avgSP);
			}
			r = rNew;
			stoppingPower_r = stoppingPower_rNew;
		}
		
	}
	
	/**
	 * Generate SEs as Poisson process over step. Any generated SEs get
	 * an equal share of the dissipated energy. The arrival rate
	 * is (1 / eps) * (electronic stopping power), where the scaling constant
	 * eps (in eV) is either given in Table 1 of
	 * (Ramachandra, Griffin, and Joy [2009]), or is estimated from it.
	 * @param r position prior to step
	 * @param rNew position after step
	 * @param stoppingPower electronic stopping power
	 */
	private void generateSEs(Coord r, Coord rNew, double stoppingPower) {
		double rate = (1.0 / EPS) * stoppingPower;
		Coord stepVector = Coord.sub(rNew, r);
		double stepLength = Coord.mag(stepVector);
		double energyToDissipate = stoppingPower * stepLength;
        double d = -log(1 - rand.nextDouble()) / rate; // "1 - " for finite d
        Electron joe;
        
        int layerIndex, numSpecies;
        double totalDensity; // in g/cm3
        double joyAvgEnergy;
		
        while (d <= stepLength) {
        	joe = new Electron();
        	joe.r = Coord.add(r, 
        			Coord.scale((d/stepLength) - 0.001, stepVector));
        	joe.rNew = Coord.add(r, Coord.scale(d/stepLength, stepVector));
        	joe.theta = PIOVERTWO;
        	joe.phi = TWOPI * rand.nextDouble();
        	
        	layerIndex = sample.getLayerIndex(joe.rNew);
            numSpecies = sample.layer[layerIndex].numSpecies;
            totalDensity = 0;
    		for (int i = 0; i < numSpecies; i++) {
    			totalDensity += 
    				sample.layer[layerIndex].species[i].speciesDensity;
    		}
    		// From Joy '95. The (* 1e3) converts from keV to eV.
        	joyAvgEnergy = pow(LAMBDA_D * totalDensity / 750, 1/1.66) * 1e3;
        	// sample from a Rayleigh distribution
        	joe.energy = ((joyAvgEnergy / sqrt(PI/2)) 
        			* sqrt(-2 * log(1-rand.nextDouble())));
        	if (energyToDissipate >= joe.energy) {
        		queue.add(joe);
        		energyToDissipate -= joe.energy;
        	} else {
        		break;
        	}
            d += -log(1 - rand.nextDouble()) / rate;
        }
        dissipateUniformly(r, rNew, energyToDissipate);
	}
	
	/**
	 * This implements the "swiss cheese" model of ion exposure,
	 * with the mesh grid dictating the width of holes in the cheese.
	 * No SEs generated. Energy is dissipated uniformly along the given
	 * leg (r to rNew) of the ion trajectory.
	 * @param r
	 * @param rNew
	 * @param stoppingPower
	 */
	private void dissipateSwissly(Coord r, Coord rNew, double stoppingPower) {
		Coord stepVector = Coord.sub(rNew, r);
		double stepLength = Coord.mag(stepVector);
		double energyToDissipate = stoppingPower * stepLength;
		dissipateUniformly(r, rNew, energyToDissipate);
	}
	

	private Particle getNextJoe() {
		return queue.poll();
	}
	
	private void calcTrajectory(Particle joe) {
		Coord rOld = joe.r;
        joe.r = joe.rNew;
		while (true) {
			/* get step length */
			int layerIndex = sample.getLayerIndex(joe.r);
			double effectiveElectronDensity = 
				sample.layer[layerIndex].effectiveElectronDensity;
			double meanExcitationEnergy =
				sample.layer[layerIndex].meanExcitationEnergy;
			int numSpecies = sample.layer[layerIndex].numSpecies;
			double[] atomDensities = new double[numSpecies];
			double[] atomicNumbers = new double[numSpecies];
			for (int i = 0; i < numSpecies; i++) {
				atomDensities[i] = 
					sample.layer[layerIndex].species[i].atomDensity;
				atomicNumbers[i] =
					sample.layer[layerIndex].species[i].atomicNumber;
			}
			double stepLength = joe.getStepLength(effectiveElectronDensity,
					atomDensities, atomicNumbers);
			/* got step length */
			
			Coord oldStepVector = Coord.sub(joe.r, rOld);
			joe.rNew = Coord.propagate(joe.r, stepLength, 
					oldStepVector, joe.theta, joe.phi);
			/*
	         * Electron passes through into next layer without scattering
	         */
	        if (joe.rNew.z > sample.layer[layerIndex].outz) {
	            double scaler = 
	            	abs((sample.layer[layerIndex].outz - joe.r.z)/
	            		(joe.rNew.z - joe.r.z));
	            joe.rNew = Coord.add(joe.r, 
	            		Coord.scale(scaler, Coord.sub(joe.rNew,joe.r)));
	            dissipate(joe,effectiveElectronDensity,meanExcitationEnergy);
	            joe.phi = 0.0;
	            joe.theta = 0.0;
	            layerIndex++;
	        }
	        /*
	         * Electron passes through into previous layer without scattering
	         */
	        else if (joe.rNew.z <= sample.layer[layerIndex].inz) {
	            double scaler = 
	            	abs((sample.layer[layerIndex].inz - joe.r.z)/
	            		(joe.rNew.z - joe.r.z));
	            joe.rNew = Coord.add(joe.r,
	            		Coord.scale(scaler, Coord.sub(joe.rNew,joe.r)));
	            dissipate(joe,effectiveElectronDensity,meanExcitationEnergy);
	            joe.phi = 0.0;
	            joe.theta = 0.0;
	            layerIndex--;
	        }
	        /*
	         * Electron scatters within current layer
	         */
	        else {
	        	dissipate(joe,effectiveElectronDensity,meanExcitationEnergy);
	        	// we are including secondaries
	        	if ((Particle.maxRank  > 0) && 
	        			// this is eligible to spawn
	        			(joe.rank != Particle.maxRank) &&
	        			joe.willSpawn(rand.nextDouble(), 
	        					effectiveElectronDensity, 
	        					atomDensities, atomicNumbers)) {
	        		/* this event is inelastic */
	        		if (Particle.cutoffE/joe.energy > 0.5) {
	        			out.println("Whoops! cutoffE is " + Particle.cutoffE
	        					+ " and particle energy is " + joe.energy);
	        		}
	        		joe = joe.spawn();
	        	} else {
	        		/* this event is elastic */
	        		if (joe.energy > 0) {
	        			joe.scatter(atomDensities, atomicNumbers);
	        		}
	        	}
	        }
	        rOld = joe.r;
	        joe.r = joe.rNew;
	        
	        if ((layerIndex < 0) || 
	        		(layerIndex >= sample.numLayers) || 
	        		(joe.energy <= 0.0)) {
	            if (joe.partner == null) {
	                if ((layerIndex < 0) && (joe.rank == 0)) {
/*	                	bseCount += 1;*/
	                }
	                break;
	            } else {
	                /*Secondary is dead, return to partner */
	            	joe = joe.partner;
	                rOld = joe.r;
	                joe.r    = joe.rNew;
	            }
	        }
			
			//break;
		}
	}

	/**
	 * Break wind
	 * @param joe
	 * @param effectiveElectronDensity
	 * @param meanExcitationEnergy
	 */
	private void dissipate(Particle joe, 
			double effectiveElectronDensity, double meanExcitationEnergy) {
		
		double deltaE,radial_r,cell_r;
	    int    i,r_index,x_index,y_index,z_index;
	    double[] Irz_ptr, Eyz_ptr, Ez_ptr;

	    Coord del   = Coord.sub(joe.rNew,joe.r);
	    double s     = Coord.mag(del);
	    int count = (int)(max(abs(del.z)/DZ, 
	    		          sqrt(del.x*del.x+del.y*del.y)/DR)) + 1;
	    del   = Coord.scale(1.0/(double)count, del);
	    double ds    = s/count;

	    Coord current = joe.r;

	    if (joe.rank != 0) {
	    	// a secondary
	        Ez_ptr = Ez_sec;
	        Irz_ptr = Irz_sec;
	        Eyz_ptr = Eyz_sec;
	    } else {
	        Ez_ptr = Ez;
	        Irz_ptr = Irz;
	        Eyz_ptr = Eyz;
	    }

	    if (joe.energy < E_MIN) {
	        /* Dump all energy in current cell */
	        radial_r = sqrt(current.x*current.x + current.y*current.y);
	        r_index = (int)(radial_r/DR);
	        x_index = (int)(abs(current.x)/DY);
	        y_index = (int)(abs(current.y)/DY);
	        z_index = (int)(current.z/DZ);
	        cell_r = ((double)r_index + 0.5)*DR;

	        Ez_ptr[z_index] += joe.energy;
	        if (z[z_index] > 0) { // we wish to track dissipation at this z
	            Irz_ptr[(z[z_index] - 1)*MAX_R_INDEX + r_index] += 
	            	joe.energy/(TWOPI*DR*DZ*cell_r);
	            Eyz_ptr[(z[z_index] - 1)*MAX_R_INDEX + x_index] += 
	            	joe.energy/(4.0*DY*DZ);
	            Eyz_ptr[(z[z_index] - 1)*MAX_R_INDEX + y_index] += 
	            	joe.energy/(4.0*DY*DZ);
	        }
	        joe.energy = 0.0;
	        return;
	    }

	    /* Step through path */
	    for (i=0;i<count;i++) {
	        deltaE  = abs(joe.dE(effectiveElectronDensity, 
	        		             meanExcitationEnergy) * ds);
	        radial_r = sqrt(current.x*current.x + current.y*current.y);
	        r_index = (int)(radial_r/DR);
	        x_index = (int)(abs(current.x)/DY);
	        y_index = (int)(abs(current.y)/DY);
	        z_index = (int)(current.z/DZ);
	        cell_r = ((double)r_index + 0.5)*DR;
	        if ((joe.energy-deltaE) < E_MIN) {
	            deltaE = joe.energy;
	            i = count; // this will be the last time through this loop
	        }
	        /* Only store energy if in region of interest */
	        if ((z_index < MAX_Z_INDEX) && (z_index >= 0)) {
	            Ez_ptr[z_index] += deltaE;
	            if ((z[z_index] > 0) && (r_index < MAX_R_INDEX)) {
	                Irz_ptr[(z[z_index] - 1)*MAX_R_INDEX + r_index] += 
	                	deltaE/(TWOPI*DR*DZ*cell_r);
	                Eyz_ptr[(z[z_index] - 1)*MAX_R_INDEX + x_index] += 
	                	deltaE/(4.0*DY*DZ);
	                Eyz_ptr[(z[z_index] - 1)*MAX_R_INDEX + y_index] += 
	                	deltaE/(4.0*DY*DZ);
	            }
	        }
	        else {
	            out.format("Out of mesh: r_index=%d, z_index=%d%n",
	            		r_index, z_index);
	        }

	        joe.energy -= deltaE;
	        current = Coord.add(current,del);

	        /* Not enough energy to get back to region of interest */
	        if (prune && (deltaE/ds*(current.z - zMax) > joe.energy/2.0)) {
	            joe.energy = 0.0;
	            i = count;
	        }
	    }
		
	}

	/**
	 * No fancy shmancy stopping power calculations here.
	 * @param r
	 * @param rNew
	 * @param energyToDissipate
	 */
	private void dissipateUniformly(Coord r, Coord rNew,
			double energyToDissipate) {
		
		// Dissipation is outside of region of interest
        if (prune && (r.z > zMax) && (rNew.z > zMax)) {return;}
		
		double deltaE,radial_r,cell_r;
	    int    i,r_index,x_index,y_index,z_index;

	    Coord del   = Coord.sub(rNew,r);
	    int count = (int)(max(abs(del.z)/DZ, 
	    		          sqrt(del.x*del.x+del.y*del.y)/DR)) + 1;
	    del   = Coord.scale(1.0/(double)count, del);

	    Coord current = r;

	    deltaE  = energyToDissipate/count;
	    /* Step through path */
	    for (i=0;i<count;i++) {
	        radial_r = sqrt(current.x*current.x + current.y*current.y);
	        r_index = (int)(radial_r/DR);
	        x_index = (int)(abs(current.x)/DY);
	        y_index = (int)(abs(current.y)/DY);
	        z_index = (int)(current.z/DZ);
	        cell_r = ((double)r_index + 0.5)*DR;
	        /* Only store energy if in region of interest */
	        if ((z_index < MAX_Z_INDEX) && (z_index >= 0)) {
	            Ez[z_index] += deltaE;
	            if ((z[z_index] > 0) && (r_index < MAX_R_INDEX)) {
	                Irz[(z[z_index] - 1)*MAX_R_INDEX + r_index] += 
	                	deltaE/(TWOPI*DR*DZ*cell_r);
	                Eyz[(z[z_index] - 1)*MAX_R_INDEX + x_index] += 
	                	deltaE/(4.0*DY*DZ);
	                Eyz[(z[z_index] - 1)*MAX_R_INDEX + y_index] += 
	                	deltaE/(4.0*DY*DZ);
	            }
	        }
	        else {
	            out.format("Out of mesh: r_index=%d, z_index=%d%n",
	            		r_index, z_index);
	        }

	        current = Coord.add(current,del);
	    }
		
	}
}
