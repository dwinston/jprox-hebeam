package sim;

public class Species {
	double atomicNumber;
	double atomicWeight;
	// density of species in layer (g/cm3)
	double speciesDensity;
	// atoms/(Angstrom)^3
	double atomDensity;
	// electrons/(Angstrom)^3
	double electronDensity;
	double excitationEnergy;
	
	/**
	 * Only update when atomicNumber, atomicWeight, 
	 * and speciesDensity are initialized.
	 */
	void update() {
		atomDensity = (Standard.N_A * speciesDensity / 
				(1.0E24 * atomicWeight));
		this.electronDensity = atomicNumber * atomDensity;
		excitationEnergy = Standard.mean_I[(int) (atomicNumber-1)];
	}
}
