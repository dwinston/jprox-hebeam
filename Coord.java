package sim;

import static java.lang.Math.sin;
import static java.lang.Math.cos;
import static java.lang.Math.sqrt;

/**
 * A 3D coordinate vector and its methods.
 * @author Donny
 *
 */
public class Coord {
	
	double x,y,z;
	
	public Coord() {
		this.x = 0;
		this.y = 0;
		this.z = 0;
	}
	
	public Coord(double x, double y, double z) {
		this.x = x;
		this.y = y;
		this.z = z;
	}
	
	static Coord add(Coord a,Coord b) {
	    Coord result = new Coord();

	    result.x = a.x + b.x;
	    result.y = a.y + b.y;
	    result.z = a.z + b.z;

	    return (result);
	}

	static Coord sub(Coord a,Coord b) {
	    Coord result = new Coord();

	    result.x = a.x - b.x;
	    result.y = a.y - b.y;
	    result.z = a.z - b.z;

	    return (result);
	}

	static Coord cross(Coord a,Coord b) {
	    Coord result = new Coord();

	    result.x = (a.y)*(b.z) - (a.z)*(b.y);
	    result.y = (a.z)*(b.x) - (a.x)*(b.z);
	    result.z = (a.x)*(b.y) - (a.y)*(b.x);

	    return (result);
	}

	static Coord scale(double a,Coord b) {
	    Coord result = new Coord();

	    result.x = a*b.x;
	    result.y = a*b.y;
	    result.z = a*b.z;

	    return (result);
	}

	static double mag(Coord a) {
	    double result;
	    result = Math.sqrt((a.x)*(a.x) + (a.y)*(a.y) + (a.z)*(a.z));
	    return (result);
	}
	
	@Override
	public boolean equals(Object o) {
		if (o instanceof Coord) {
			Coord c = (Coord) o;
			if (c.x == x && c.y == y && c.z == z) {
				return true;
			}
		}
		return false;		
	}
	
	@Override
	public int hashCode() {
		return (int) (100000*x + 1000*y + z); 
	}

	public static Coord propagate(Coord r, double stepLength,
			Coord oldStepVector, double theta, double phi) {
		
		Coord result;
		double sintcosp, sintsinp, cost, root;
		Coord l = oldStepVector;
		Coord l_prime = new Coord(); // "newStepVector"
		
		sintcosp = sin(theta)*cos(phi);
	    sintsinp = sin(theta)*sin(phi);
	    cost = cos(theta);

	    l = scale(1.0/mag(l),l);

	    root = sqrt(l.x*l.x + l.y*l.y);

	    if (root > 0.0)
	    {
	        l_prime.x = l.x*cost + l.y/root*sintcosp + l.x*l.z/root*sintsinp;
	        l_prime.y = l.y*cost - l.x/root*sintcosp + l.y*l.z/root*sintsinp;
	        l_prime.z = l.z*cost - root*sintsinp;
	    }
	    else
	    {
	        l_prime.x = sintcosp;
	        l_prime.y = sintsinp;
	        l_prime.z = cost;
	    }

	    l_prime = scale(stepLength,l_prime);
	    result = add(r,l_prime);

	    return (result);
	}

}
