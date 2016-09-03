/*
 * Based on the following paper:
 * Geoff Hulten, Laurie Spencer and P. Domingos
 * "Mining  Time-Changing Data  Streams",  
 * In  Proceedings  of  the  ACM  SIGKDD 2001. 
 * 
 * Description:
 * Each concept is modeled by alternating class bands separated by parallel hyper-planes.
 * Details in section 4 of the paper.
 * 
 * @author: Ariyam Das
 */


import java.util.Random;

public class SyntheticConceptForDrift {
	
	final private int NUM_CLASSES = 2;
	
	private int dims;		//number of attributes
	private double wts[];	//weights indicating information content of different attributes
	private double w0;		//hyper-plane constant
	
	//initialization done according to the details in the paper
	public SyntheticConceptForDrift(int dims)
	{
		this.dims = dims;
		
		wts = new double[dims];
		for(int i=0; i<dims; i++)
			wts[i] = 0.2;
		
		w0 = 0.25*dims;	
	}
	
	public int getDims()
	{
		return this.dims;
	}
	
	public double [] getHyperplaneWts()
	{
		return this.wts;
	}
	
	public double getHyperplaneWt(int i)
	{
		return this.wts[i];
	}
	
	public double getHyperplaneConst()
	{
		return this.w0;
	}
	
	//returns the class label for the given input point based on the current concept.
	public int classify(double pt[])
	{
		if(pt.length < dims)
		{
			System.out.println("Dimensionality of the point less than hyperplane!");
			System.exit(-1);
		}
		
		double s=0;
		for(int i=0; i<dims; i++)
			s += wts[i]*pt[i];
		
		int band = (int)(Math.abs(s)/w0 * 10) % NUM_CLASSES;
		return band;
	}
	
	//returns a random class label
	public int randomClassify()
	{
		Random rand = new Random();
		return rand.nextInt(NUM_CLASSES);
	}
	
	public void updateHyperplaneWt(int i, double delta)
	{
		this.wts[i] += delta;
	}
	
	public void updateHyperplaneWts(double delta)
	{
		for(int i=0; i<this.dims; i++)
			this.wts[i] += delta;
	}
	
	public void updateHyperplaneWts(int index[], double delta[])
	{
		for(int i=0; i<index.length; i++)
			this.wts[index[i]] += delta[i];
	}
	
	@Override
	public String toString()
	{
		String string = "";
		string = "Dims:" + dims + "\n" + "w0:" + w0 + "\n" + "Wts:";
		for(int i=0; i<dims; i++)
			string += wts[i]+" ";
			
		return string;
	}
	
	public static void main(String args[]) 
	{
		SyntheticConceptForDrift scd = new SyntheticConceptForDrift(5);
		System.out.println(scd);	
	}
}
