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

import java.util.Arrays;
import java.util.Hashtable;
import java.util.ArrayList;
import java.util.List;
import java.util.Random;
import java.io.File;
import java.io.FileInputStream;
import java.io.BufferedReader;
import java.io.InputStreamReader;
import java.io.IOException;
import java.util.stream.*;

public class HoeffdingTreeBLB {
	
	//following should have been pushed into the constructor;
	//but for sake of simplicity and all practical implementation purposes, this will do!
	//private final int MAX_CLASS = 2;	
	private final int MAX_ATT_VAL = 2;
	
	private final double USER_DEFINED_THRESHOLD = 0.1;
	private final int REBOOT_CUT = 100;
	
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
		if(ht.runBLB() && !ht.isPure())
		{	
			//run BLB - ARG : all_pts arraylist
			//what is the estimator vector? -- call the method entropyVector for each sample. 
			//each sample is also an array list
			//from each sample : Vector has N+1 estimators where N is the number of dimensions
			//average all vector results
			//avg_ent_vector
			double[] avg_ent_vector;
			//Test, exactly same as HoeffdingTree (without BLB)
			//if (Math.pow(data.length, .5) * 3 * 30 >= data.length) {
	//		avg_ent_vector = ht.entropyVector(ht.all_pts);
	//			System.out.println("b");
			//}
			//======================================
			//else {
				//System.out.println(ht.all_pts.size());
				//System.out.println(ht.all_pts.size());
				//System.out.println((int) Math.pow(data.length, .8));
			int s = 4;
			int b = (int) Math.pow(data.length,.6);
			int r = 75;
			int n = ht.all_pts.size();
			double[] jmst = new double[ht.all_pts.get(0).length];
			for (int j=0;j<s;j++) {
				int[] partition = fyShuffle(s,b,n);
				double[] kmst = new double[ht.all_pts.get(0).length];
				for (int k=0;k<r;k++) {	
					kmst = sumVector(ht.entropyVector(multiDist(ht.all_pts,partition,s,b)),kmst);
				}
				jmst = sumVector(avgVector(kmst,r),jmst);
			}
			avg_ent_vector = avgVector(jmst,s);
				//avg_ent_vector = ht.runBLB(ht.all_pts, 4, (int) Math.pow(data.length,.5) , 40);
				//System.out.println("a");
			//}
			//bestSplittingAttr(avg_ent_vector)
			//whether to split or not
			int best_attr = ht.getBestSplitAttribute(avg_ent_vector,ht.pts_read);
			if(best_attr != -1) 
				ht.split(best_attr); 
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
            			System.out.print(no_of_eg+","+(double)corr/te+","+this.statCountNodes()+";");
            			
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
	
	public void split(int by_attr)
	{
		attribute = by_attr;
		children_cnt = MAX_ATT_VAL;
		edge_label = new int[children_cnt];
		children = new HoeffdingTreeBLB[children_cnt];
		
		for(int i=0; i<children_cnt; i++)
		{
			edge_label[i] = i;
			
			boolean pool[] = attr_set.clone();
			pool[by_attr] = true;
			
			int cl = getDomClass(by_attr,i);
			cl = ((cl != -1) ? cl : majority_class);
			
			children[i] = new HoeffdingTreeBLB(delta, pool, cl);
		}	
		//clear up space
		attr_set = null;	//gc should clear this
		all_pts.clear();	//cleared
		all_pts = null;		//gc will clear this now
	}
	
	//will work for binary classes only for now
	public int getDomClass(int id, int val)
	{
		int ncl = 0, ycl = 0;
		for(int i=0; i<all_pts.size(); i++)
		{
			String data[] = all_pts.get(i);
			if(Integer.parseInt(data[id])==val)
			{
				if(Integer.parseInt(data[data.length-1])==0)
					ncl++;
				else
					ycl++;
			}
		}
		return ((ncl>ycl)?0:1);
	}

	//BLB Section
	//TODO replace with proper distribution, temp naive algorithm for speed
	private int[] fyShuffle(int s, int b, int row) {
		boolean[] check = new boolean[row];
		int[] res = new int[b];
		Random r = new Random();
		int cnt = 0;
		while(cnt<b) {
			int ind = r.nextInt(row);
			if (check[ind]==false) {
				check[ind] = true;
				res[cnt] = ind;
				cnt++;
			}
		}
		return res;
	}
	
	private ArrayList<String[]> multiDist(ArrayList<String[]> dataset, int[] partition, int s, int len) {
		int n = dataset.size();
		ArrayList<String[]> ret = new ArrayList<String[]>();

		Random rand = new Random();
		int remain = n;
		while (remain > 0) {
			int index = rand.nextInt(len);
			int amt = rand.nextInt(remain)+1;
			for (int i=0;i<amt;i++) {
				ret.add(dataset.get(partition[index]));
			}
			remain -= amt;
		}
		return ret;
	}
	
	private double[] sumVector(double[] in, double[] old) {
		double[] sum = new double[old.length];
		if(in.length != old.length) {
			return old;
		}
		else {
			for(int i=0;i<old.length;i++) {
				sum[i] = in[i] + old[i];
			}
			return sum;
		}		
	}
	
	private double[] avgVector(double[] in, int num) {
		double[] avg = in;
		for (int i=0;i<in.length;i++) {
			avg[i] = in[i]/num;
		}
		return avg;
	}
	
	
	
	public static void main(String args[]) 
	{
		HoeffdingTreeBLB ht = new HoeffdingTreeBLB(Math.pow(10,-7), 100);	//delta=0.1%, 100 attributes
		//ht.streamInputFile("data/0.25_0.2_10241_5121.dat");
		
		ht.streamInputFile("",10000,1000);
		//System.out.println(ht.statCountNodes());
		//System.out.println(ht.min_pts);
		//for(int i=0; i<100; i++)
		//System.out.println(i+"-"+ht.attributeEntropyCal(i, ht.pts_read));
		//System.out.println(ht.getBestSplitAttribute());
		//System.out.println(ht.statCountNodes());
		//System.out.println(ht.statCountLeaves());
	}
}
