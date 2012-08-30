package sim;

public class Layer {
	double thickness;
	int numSpecies;
	Species[] species;
	double meanExcitationEnergy;
	double effectiveElectronDensity;
	double inz, outz;
	
	/**
	 * Only update when species array is initialized and set.
	 * Not required for update: thickness, inz, outz.
	 */
	void update() {
		numSpecies = species.length;
		for (int i = 0; i < numSpecies; i++) {
			effectiveElectronDensity += species[i].electronDensity;
			meanExcitationEnergy += (species[i].electronDensity * 
					Math.log(species[i].excitationEnergy));
		}
		meanExcitationEnergy /= effectiveElectronDensity;
		meanExcitationEnergy = Math.exp(meanExcitationEnergy);
	}
}
