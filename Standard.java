package sim;

/**
 * Provides a few static final methods and constants for general use.
 * @author Donny
 *
 */
public class Standard {
	public static final double PI = Math.PI;
	public static final double TWOPI = 2 * PI;
	public static final double PIOVERTWO = PI / 2;
	public static final double Q = 1.60219E-19;
	public static final double ESU = 4.80324E-10;
	public static final double N_A = 6.022E23;
	/* scaling constant for SE generation rate of HeIons, 
	 * from Ramachandra 2009. */
	public static final double EPS = 60.0; // in eV
	/* escape depth for SEs, from Ramachandra 2009. 
	 * Varies with atomic number Z, but there are no tabulated values 
	 * for Si, H, or O, nor for any compounds, so here I make a guess 
	 * for both Si and HSQ: the minimum plus the half-range for Z > 3. */
	public static final double LAMBDA_D = 9.75; // in Angstroms
	
	/**
	 * Table of mean ionization energies for the elements.
	 * Source: NIST
	 * http://physics.nist.gov/PhysRefData/XrayMassCoef/tab1.html
	 * JF 11/15/08
	 */
	public static final double[] mean_I = 
	{19.2, 41.8, 40, 63.7, 76, 78, 82, 95,
	115, 137, 149, 156, 166, 173, 173, 180, 174, 188, 190, 191, 216, 233,
	245, 257, 272, 286, 297, 311, 322, 330, 334, 350, 347, 348, 343, 352,
	363, 366, 379, 393, 417, 424, 428, 441, 449, 470, 470, 469, 488, 488,
	487, 485, 491, 482, 488, 491, 501, 523, 535, 546, 560, 574, 580, 591,
	614, 628, 650, 658, 674, 684, 694, 705, 718, 727, 736, 746, 757, 790,
	790, 800, 810, 823, 823, 830, 825, 794, 827, 826, 841, 847, 878, 890 };
	
}
