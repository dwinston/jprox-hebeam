package sim;

public class Sample {
	int numLayers;
	Layer[] layer;
	
	public Sample(int numLayers) {
		this.numLayers = numLayers;
		layer = new Layer[numLayers];
	}
	
	/**
	 * Only update when layer array is initialized and 
	 * each layer has set thickness.
	 */
	public void update() {
		double z = 0;
		for (Layer L : layer) {
			L.inz = z;
			L.outz = L.inz + L.thickness;
			z = L.outz;
		}
	}
	
	int getLayerIndex(Coord r) {
		double z = r.z;
		int index = -1;
		for (Layer L : layer) {
			if (z < L.inz) {return index;} 
			else {index += 1;}
		}
		if (z < layer[layer.length - 1].outz) {return index;} 
		else {return index+1;}
	}
}
