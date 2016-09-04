/*
 * Based on the following paper:
 * Geoff Hulten, Laurie Spencer and P. Domingos
 * "Mining  Time-Changing Data  Streams",  
 * In  Proceedings  of  the  ACM  SIGKDD 2001.
 * 
 * Description:
 * Synthetic data generation in presence of concept drift.
 * See details from the paper. 
 * 
 * @author: Ariyam Das
 */

import java.util.ArrayList;
import java.util.Random;
import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.IOException;

public class SyntheticConceptDriftDataGeneration {
	
	private int dims;
	private int total_pts;
	private int drift_at;
	private int drifting_attributes;
	private double noise;
	private SyntheticConceptForDrift target;
	private String output_filename;
	
	public SyntheticConceptDriftDataGeneration(int d, int num_pts, int drft_at, double n)
	{
		this.dims = d;
		target = new SyntheticConceptForDrift(dims);	
		this.total_pts = num_pts;
		this.drift_at = drft_at;
		this.drifting_attributes = 1;	//default value
		this.noise = n;	
		output_filename = "data/d-"+dims+"_n-"+noise+"_pts-"+total_pts+"_at-"+drift_at+"_dr-"+drifting_attributes+".dat";		
		generateDataFile();				//stored according to the output filename
	}
	
	public SyntheticConceptDriftDataGeneration(int d, int num_pts, int drft_at, int drft_level, double n)
	{
		this.dims = d;
		target = new SyntheticConceptForDrift(dims);
		this.total_pts = num_pts;
		this.drift_at = drft_at;
		this.drifting_attributes = drft_level;	
		this.noise = n;		
		output_filename = "data/d-"+dims+"_n-"+noise+"_pts-"+total_pts+"_at-"+drift_at+"_dr-"+drifting_attributes+".dat";		
		generateDataFile(drifting_attributes);		
	}
		
	public String generateRandomTrainingInstance()
	{
		Random rand = new Random();
		double point[] = new double[dims+1];
		
		for(int i=0; i<dims; i++)	
			point[i] = rand.nextDouble();
		
		if(Math.random() < noise)
			point[dims] = target.randomClassify(); //class affected by noise
		else
			point[dims] = target.classify(point);	//assign class
		
		String string = ""+point[0];
		for(int i=1; i<point.length; i++)
		{
			string += "," + point[i];
		}
		return string;	
	}
	
	//Concept drift based data generated based on one hyper-plane parameter, w1.
	public void generateDataFile()
	{
		try 
    	{
			BufferedWriter out = new BufferedWriter(new FileWriter(output_filename,false));
			
			//No concept drift for the first set .....
			for(int x=0; x<drift_at; x++)
			{
				String text = generateRandomTrainingInstance();
	    		out.append(text+"\n");	
			}
			
			//Now the concept drift begins .....
			double delta = 0.01*dims*1;
			int flip_at=5;
			for(int x=drift_at; x<total_pts; x+=drift_at)
			{
				double val = target.getHyperplaneWt(0);

				if((x/drift_at)%flip_at == 0)
					delta *= -1;			
				
				if(val+delta<0 || val+delta > 0.25*dims)
					delta *= -1;
				
				target.updateHyperplaneWt(0, delta);		//adding concept drift
				System.out.println("Drifted w1: "+ target.getHyperplaneWt(0));
				
				for(int y=0; y<drift_at; y++)
				{
					String text = generateRandomTrainingInstance();
		    		out.append(text+"\n");	
				}	
			}
    		out.close();
    	}
    	catch ( IOException e )
    	{
    		e.printStackTrace();
    	}
	}
	
	//Concept drift data generated based on D randomly chosen hyper-plane parameters.
	public void generateDataFile(int D)
	{
		try 
    	{
			BufferedWriter out = new BufferedWriter(new FileWriter(output_filename,false));
			
			//No concept drift for the first set .....
			for(int x=0; x<drift_at; x++)
			{
				String text = generateRandomTrainingInstance();
	    		out.append(text+"\n");	
			}
			
			//Now the concept drift begins .....
			double delta [] = new double[dims]; 
			for(int i=0; i<dims; i++)
				delta[i] = 0.01*dims*1;
			
			for(int x=drift_at; x<total_pts; x+=drift_at)
			{
				int param[] = randomGenerate(D);	//randomly chosen D hyper-plane parameters
				for(int i=0; i<D; i++)
				{
					double val = target.getHyperplaneWt(param[i]);

					Random rand = new Random();
					if(rand.nextInt(4)==0)
						delta[param[i]] *= -1;	//25% chance of sign flip
					
					if(val+delta[param[i]] < 0 || val+delta[param[i]] > 0.25*dims)
						delta[param[i]] *= -1;
					
					target.updateHyperplaneWt(param[i], delta[param[i]]);	//adding concept drift	
					System.out.println("Drifted w"+param[i]+": "+ target.getHyperplaneWt(param[i]));
				}
				
				for(int y=0; y<drift_at; y++)
				{
					String text = generateRandomTrainingInstance();
		    		out.append(text+"\n");	
				}	
			}
    		out.close();
    	}
    	catch ( IOException e )
    	{
    		e.printStackTrace();
    	}
	}
	
	public int [] randomGenerate(int D)
	{
		int vals [] = new int[D];
		int nums [] = new int[dims];
		for(int i=0; i<dims; i++)
			nums[i]=i;
		Random rand = new Random();
		for(int i=0; i<D; i++)
		{
			int x = rand.nextInt(dims-i);
			vals[i] = nums[x];	//the random value
			nums[x] = nums[dims-i-1];	//bring in the last value
		}
		return vals;
	}
	
	public static void main(String args[]) 
	{
		//<#dims, #total points, #drift at every this many points, #noise level>
		new SyntheticConceptDriftDataGeneration(25, 500000, 5000, 0.1);
		
		//<#dims, #total points, #drifting points, #drifting attributes, #noise level>
		//new SyntheticConceptDriftDataGeneration(25, 500000, 5000, 5, 0.1);
	}
}
