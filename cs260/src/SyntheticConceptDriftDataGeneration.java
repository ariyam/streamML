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
		generateDataFile(drifting_attributes);		//stored according to the output filename
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
	
	public void generateDataFile()
	{
		try 
    	{
			BufferedWriter out = new BufferedWriter(new FileWriter(output_filename,true));
			
			//No concept drift for the first set .....
			for(int x=0; x<drift_at; x++)
			{
				String text = generateRandomTrainingInstance();
	    		out.append(text+"\n");	
			}
			
			//Now the concept drift begins .....
			for(int x=drift_at; x<total_pts; x+=drift_at)
			{
				//concept drift - code
				double delta = 0.01*dims*1;			//from section 4 in paper
				
				if ((x/drift_at)%2 == 1)		//keep alternating between two weights for w1
					target.updateHyperplaneWt(0, delta);	//updated the weight
				else
					target.updateHyperplaneWt(0, -delta);	//revert back
				
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
			for(int x=drift_at; x<total_pts; x+=drift_at)
			{
				//concept drift - code
				double delta = 0.01*dims*1;			//from section 4 in paper
				
				if ((x/drift_at)%2 == 1)		//keep alternating between two weights for w1
					{target.updateHyperplaneWt(0, delta);	//updated the weight
					System.out.println(target.getHyperplaneWt(0));}
				else
					{target.updateHyperplaneWt(0, -delta);	//revert back
					System.out.println(target.getHyperplaneWt(0));}
					
				
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
	
	
	public static void main(String args[]) 
	{
		//<#dims, #total points, #drift at every this many points, #noise level>
		new SyntheticConceptDriftDataGeneration(25, 500000, 5000, 1, 0.1);
	}
}
