/*
 * Based on the following paper:
 * P. Domingos and G. Hulten,  
 * "Mining  high-speed  data  streams",  
 * In  Proceedings  of  the  ACM  SIGKDD 2000. 
 * 
 * Description:
 * Hoeffding tree generated from the sequence of incoming examples,
 * following the algorithm mentioned in the paper.
 * 
 * @author: Ariyam Das
 * 
 * Note: n_min for leaves (minimum number of points before G is computed),
 * as mentioned in the paper, has not been implemented yet.
 * 
 */

import java.util.Hashtable;
import java.util.ArrayList;
import java.io.File;
import java.io.FileInputStream;
import java.io.BufferedReader;
import java.io.InputStreamReader;
import java.io.IOException;

public class HoeffdingTreeBLB {
	
	//following should have been pushed into the constructor;
	//but for sake of simplicity and all practical implementation purposes, this will do!
	private final int MAX_CLASS = 2;	
	private final int MAX_ATT_VAL = 2;
	
	private final double USER_DEFINED_THRESHOLD = 0.05;
	private final int REBOOT_CUT = 1;
	
	private int attribute;
	private int children_cnt;
	private int majority_class;
	private int majority_votes;
	private int [] edge_label;
	private HoeffdingTreeBLB [] children;
	private double delta;
	private double R;
	private int pts_read;
	boolean []attr_set;
	ArrayList<String[]> all_pts;
	private int reboot_pts;
	
	public HoeffdingTreeBLB(double delta, int no_attr)
	{
		attr_set = new boolean[no_attr];	//initialized to false by default
		this.delta = delta;
		this.R = 1;
		this.children_cnt = 0;
		this.attribute = -1;
		this.majority_class = -1;
		this.majority_votes = 0;
		this.pts_read = 0;
		this.edge_label = null;
		this.children = null;
		this.all_pts = new ArrayList<String[]>();
		this.reboot_pts = 0;
	}
	
	public HoeffdingTreeBLB(double delta, double R, int no_attr)
	{	
		attr_set = new boolean[no_attr];	//initialized to false by default
		this.delta = delta;
		this.R = R;
		this.children_cnt = 0;
		this.attribute = -1;
		this.majority_class = -1;
		this.majority_votes = 0;
		this.pts_read = 0;
		this.edge_label = null;
		this.children = null;
		this.all_pts = new ArrayList<String[]>();
		this.reboot_pts = 0;
	}
	
	public HoeffdingTreeBLB(double delta, boolean unsel_pool[], int maj_class)
	{
		this.delta = delta;
		this.R = 1;
		this.children_cnt = 0;
		this.attribute = -1;
		this.majority_class = maj_class;
		this.majority_votes = 0;
		this.pts_read = 0;
		this.edge_label = null;
		this.children = null;
		this.attr_set = unsel_pool;
		this.all_pts = new ArrayList<String[]>();
		this.reboot_pts = 0;
	}
	
	public HoeffdingTreeBLB(double delta, double R, boolean unsel_pool[], int maj_class)
	{
		this.delta = delta;
		this.R = R;
		this.children_cnt = 0;
		this.attribute = -1;
		this.majority_class = maj_class;
		this.majority_votes = 0;
		this.pts_read = 0;
		this.edge_label = null;
		this.children = null;
		this.attr_set = unsel_pool;
		this.all_pts = new ArrayList<String[]>();
		this.reboot_pts = 0;
	}
	
	public boolean isLeaf()
	{
		if(this.children_cnt == 0)
			return true;
		else
			return false;
	}
	
	public int getAttribute()
	{
		return this.attribute;
	}
	
	public int predictClass()
	{
		return this.majority_class;
	}
	
	public int predictClass(String tuple[])
	{
		HoeffdingTreeBLB ht = this;
		while(!ht.isLeaf())
		{
			int a_no = ht.getAttribute();
			int e_val = Integer.parseInt(tuple[a_no]);
			ht = ht.getChild(e_val);
		}

		return ht.predictClass();
	}
	
	public void updateClass(int cl)
	{
		if(majority_votes>0) 
		{
			if(cl == majority_class)
				majority_votes++;
			else
				majority_votes--;
		}
		else
		{
			majority_class = cl;
			majority_votes++;
		}
	}
	
	public HoeffdingTreeBLB getChild(int edge_val)
	{
		for(int i=0; i<children_cnt; i++)
			if(edge_label[i]==edge_val)
				return children[i];
		return null;
	}
	
	public boolean isPure()
	{
		return (majority_votes == pts_read);
	}
	
	public boolean runBLB()
	{
		if (reboot_pts == REBOOT_CUT)
		{
			reboot_pts = 0;
			return true;
		}
		return false;
	}
	
	public void updateParams(String data[])
	{
		pts_read++;
		
		int cols = data.length - 1;
		int cl = Integer.parseInt(data[cols]);
		updateClass(cl);
		
		reboot_pts++;
		
		all_pts.add(data);		
	}
	
	public void processTuple(String data[])
	{
		HoeffdingTreeBLB ht = this;
		while(!ht.isLeaf())
		{
			int dim_no = ht.getAttribute();
			ht = ht.getChild(Integer.parseInt(data[dim_no]));
		}
		
		ht.updateParams(data);
		if(!ht.isPure() && ht.runBLB())
		{
			//run BLB
			//what is the estimator?
		}
	}
	
	public void streamInputFile(String filename)
	{
		try 	//read input file
		{
			File file = new File(filename);
			FileInputStream fis = new FileInputStream(file);
			BufferedReader br = new BufferedReader(new InputStreamReader(fis));
    		String line = null;
            
            while ((line = br.readLine()) != null) 
            {
            	String data[] = line.split(",");
            	processTuple(data);
        	}	
            br.close();          
    	}
    	catch(IOException ioe)
    	{
    		ioe.printStackTrace();
    	}	
	}
	
	public void streamInputFile(String filename, int train_size, int test_size)
	{
		try 	//read input file
		{
			File file = new File(filename);
			FileInputStream fis = new FileInputStream(file);
			BufferedReader br = new BufferedReader(new InputStreamReader(fis));
    		String line = null;
    		
            int tr = 0, te = 0, corr = 0, no_of_eg=0;          
            while ((line = br.readLine()) != null) 
            {
            	String data[] = line.split(",");
            	if(tr<train_size)
            	{
            		tr++;
            		processTuple(data);
            	}
            	else
            	{
            		te++;
            		int orig_cl = Integer.parseInt(data[data.length-1]);
            		int pred_cl = predictClass(data);
            		if(orig_cl == pred_cl)
            			corr++;
            		
            		if(te==test_size)
            		{
            			no_of_eg+=tr;
            			System.out.println(no_of_eg+","+(double)corr/te);
            			
            			tr=0;	//all set back to zero
            			te=0;
            			corr=0;
            		}
            	}
        	}	
            br.close();          
    	}
    	catch(IOException ioe)
    	{
    		ioe.printStackTrace();
    	}	
	}
	
	public long statCountNodes()
	{
		if(children_cnt==0)
			return 1;
		
		int res = 1;
		for(int i=0; i<children_cnt; i++)
			res += children[i].statCountNodes();
		return res;	
	}
	
	public long statCountLeaves()
	{
		if(children_cnt==0)
			return 1;
		
		int res = 0;
		for(int i=0; i<children_cnt; i++)
			res += children[i].statCountLeaves();
		return res;	
	}
	
	public double getEpsilon(int n)
	{
		if(pts_read==0)
			return 999999.0;
		
		double term = -R*R*Math.log(delta)/(2*n);
		double eps = Math.sqrt(term);
		return eps;
	}
	
	public double[] entropyVector(ArrayList<String[]> dataset)
	{
		int CODE_CLS_CNT = -1;
		int CODE_ATT_VAL_CNT = -2;
		int CODE_ATT_VAL_CLS_CNT = -3;
		
		Hashtable<String,Integer> params = new Hashtable<String,Integer>();
		
		for(int n=0; n<dataset.size(); n++)
		{
			String data[] = dataset.get(n);
			
			int cols = data.length - 1;
			int cl = Integer.parseInt(data[cols]);
			String cl_code = CODE_CLS_CNT+":"+cl;
			if(params.containsKey(cl_code))
			{
				int x = params.get(cl_code);
				params.put(cl_code, x+1);
			}
			else
				params.put(cl_code, 1);

			for(int i=0; i<cols; i++)
			{
				int val = Integer.parseInt(data[i]);
				String at_code = CODE_ATT_VAL_CNT+":"+i+":"+val;
				String at_cl_code = CODE_ATT_VAL_CLS_CNT+":"+i+":"+val+":"+cl;
				
				if(params.containsKey(at_code))
				{
					int x = params.get(at_code);
					params.put(at_code, x+1);
				}
				else
					params.put(at_code, 1);
				
				if(params.containsKey(at_cl_code))
				{
					int x = params.get(at_cl_code);
					params.put(at_cl_code, x+1);
				}
				else
					params.put(at_cl_code, 1);
			}	//end of reading tuple	
		}	//end of reading bootstrap sample
		
		double ret_vector[] = new double[dataset.get(0).length];
		for(int i=0; i<ret_vector.length-1; i++)
			if(!attr_set[i])
				ret_vector[i] = attributeEntropyCal(params, i,dataset.size());
			else
				ret_vector[i] = 100;
		
		String cl_code = CODE_CLS_CNT+":"+0;
		int n_cl = ((params.containsKey(cl_code)) ? params.get(cl_code) : 0);
		int y_cl = dataset.size() - n_cl;
		ret_vector[ret_vector.length-1] = binaryEntropyCal(y_cl,n_cl);
		
		return ret_vector;
	}
	
	public double binaryEntropyCal(int x, int y)
	{
		double res = 0;
		double s = x+y;
		
		if(x==0 || y==0)
			return res;
		
		res = -((double)x/s)*(Math.log(x/s)/Math.log(2)) -((double)y/s)*(Math.log(y/s)/Math.log(2));
		//res = 1 - (x/s)*(x/s) - (y/s)*(y/s);	//Gini index
		return res;
	}
	
	//For now, will work for binary classes and binary attributes only
	public double attributeEntropyCal(Hashtable<String, Integer> params, int att_id, int cnt)
	{
		int CODE_ATT_VAL_CLS_CNT = -3;
		
		double res = 0;
		
		for(int att_val=0; att_val<MAX_ATT_VAL; att_val++)
		{
			String code = CODE_ATT_VAL_CLS_CNT+":"+att_id+":"+att_val+":"+0;
			int n_cl = ((params.containsKey(code)) ? params.get(code) : 0);
				
			code = CODE_ATT_VAL_CLS_CNT+":"+att_id+":"+att_val+":"+1;
			int y_cl = ((params.containsKey(code)) ? params.get(code) : 0);
				
			double wt = (n_cl+y_cl);
			res += wt * binaryEntropyCal(n_cl, y_cl);
		}
		res/=cnt;
		return res;
	}
	
	//For now, will work for binary classes and binary attributes only
	public int getBestSplitAttribute(double[] ent_vector, int nval)
	{
		int b_attr = -1;
			
		//highest inf. gain means lowest entropy
		double fmin_ent = 100; //1st minimum entropy
		double smin_ent = 100; //2nd minimum entropy
		int tmp=-1;
			
		double nsplit_ent = ent_vector[ent_vector.length-1]; //entropy for no split

		for(int i=0; i<attr_set.length; i++)
		{
			if(!attr_set[i])	//this attribute not used in path from root to this node
			{
				double x = ent_vector[i];
				if(x<fmin_ent)
				{
					smin_ent = fmin_ent;
					fmin_ent = x;
					tmp = i;
				}
				else
				{
					if(x<smin_ent)
						smin_ent = x;
				}
			}
		}
			
		double x = Math.min(smin_ent, nsplit_ent);
		//double diff = (x-fmin_ent)/(nsplit_ent - fmin_ent); //percentage points did not work
		double diff = (x-fmin_ent); //simple difference, according to algo in paper
		double eps = getEpsilon(nval);
			
			if(diff > eps || (diff < eps && USER_DEFINED_THRESHOLD > eps))	//avoid comparing two very close attributes
				b_attr = tmp;
			
			return b_attr;
		}

	//BLB Section
	/*
	private static int[] fyShuffle(int[] init) {
		Random r = new Random();
		int maxbound = init.length;
		for (int j=maxbound-1;j>0;j--) {
			int k = r.nextInt(j+1);
			int temp = init[k];
			init[k] = init[j];
			init[j] = temp;
		}
		return init;
	}
	
	private static int[] multiDist(int len) {
		int ret[] = new int[len];
		Random rand = new Random();
		int remain = len;
		int ind = 0;
		while(remain > 0) {
			int j = rand.nextInt(remain);
			remain-=j;
			remain+=ret[ind%len];
			ret[ind%len] = j;
			ind++;
		}
		ret = fyShuffle(ret);
		return ret;
	}
	
	//TODO Threaded runs
	public boolean runBLB(String[] data, int s, int b, int r) {
		HoeffdingTree ht = this;
		int n = data.length;
		int[] init = new int[n];
		for (int in=0;in<n;in++) {
			init[in] = in;
		}
		int[] partition = fyShuffle(init);
		
		for (int j=0;j<s;j++) {
			int qsum = 0;
			for (int k=0;k<r;k++) {
				int[] dist = multiDist(b);
				for (int p=0;p<b;p++) {
					//TODO write function to input data
					int best_attr = ht.getBestSplitAttribute();
				}
			}
		}
		
		return true; //TODO incomplete. In progress
	}
	*/
	
	public static void main(String args[]) 
	{
		HoeffdingTreeBLB ht = new HoeffdingTreeBLB(Math.pow(10,-7), 100);	//delta=0.1%, 100 attributes
		//ht.streamInputFile("data/0.25_0.2_10241_5121.dat");
		ht.streamInputFile("data/0.15_0.0_70127_35064.dat",100000,10000);
		//System.out.println(ht.statCountNodes());
		//System.out.println(ht.min_pts);
		//for(int i=0; i<100; i++)
		//System.out.println(i+"-"+ht.attributeEntropyCal(i, ht.pts_read));
		//System.out.println(ht.getBestSplitAttribute());
		//System.out.println(ht.statCountNodes());
		//System.out.println(ht.statCountLeaves());
	}
}
