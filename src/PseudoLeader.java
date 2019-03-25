import comp34120.ex2.PlayerImpl;
import comp34120.ex2.PlayerType;
import comp34120.ex2.Record;
import java.rmi.NotBoundException;
import java.rmi.RemoteException;
import java.util.List;
import java.util.ArrayList;
import java.lang.Math;

/**
 * A pseudo leader. The members m_platformStub and m_type are declared
 * in the PlayerImpl, and feel free to use them. You may want to check
 * the implementation of the PlayerImpl. You will use m_platformStub to access
 * the platform by calling the remote method provided by it.
 * @author Xin
 */
final class PseudoLeader
	extends PlayerImpl
{

	List<double[][]> oldParams = new ArrayList<double[][]>();
	List<double[][]> oldP = new ArrayList<double[][]>();
	static int currentDate;
	static int steps;


	final double LAMBDA = 0.95;
	final double COST = 1.0;

	

	/**
	 * In the constructor, you need to call the constructor
	 * of PlayerImpl in the first line, so that you don't need to
	 * care about how to connect to the platform. You may want to throw
	 * the two exceptions declared in the prototype, or you may handle it
	 * by using "try {} catch {}". It's all up to you.
	 * @throws RemoteException
	 * @throws NotBoundException
	 */
	PseudoLeader()
		throws RemoteException, NotBoundException
	{
		/* The first parameter *MUST* be PlayerType.LEADER, you can change
		 * the second parameter, the name of the leader, such as "My Leader" */
		super(PlayerType.LEADER, "Pseudo Leader");
		//oldParams = new ArrayList<double[][]>();
		
	}

	/**
	 * You may want to delete this method if you don't want modify
	 * the original connection checking behavior of the platform.
	 * Actually I recommend you to delete this method from your own code
	 * @throws RemoteException If implemented, the RemoteException *MUST* be
	 * thrown by this method
	 */
	@Override
	public void checkConnection()
		throws RemoteException
	{
		super.checkConnection();
		//TO DO: delete the line above and put your own code here
	}

	/**
	 * You may want to delete this method if you don't want the platform
	 * to control the exit behavior of your leader class
	 * @throws RemoteException If implemented, the RemoteException *MUST* be
	 * thrown by this method
	 */
	@Override
	public void goodbye()
		throws RemoteException
	{
		super.goodbye();
		//TO DO: delete the line above and put your own exit code here
	}

	/**
	 * You may want to delete this method if you don't want to do any
	 * initialization
	 * @param p_steps Indicates how many steps will the simulation perform
	 * @throws RemoteException If implemented, the RemoteException *MUST* be
	 * thrown by this method
	 */
	@Override
	public void startSimulation(int p_steps)
		throws RemoteException
	{
		steps = p_steps;
		super.startSimulation(p_steps);
		//TO DO: delete the line above and put your own initialization code here
	}

	/**
	 * You may want to delete this method if you don't want to do any
	 * finalization
	 * @throws RemoteException If implemented, the RemoteException *MUST* be
	 * thrown by this method
	 */
	@Override
	public void endSimulation()
		throws RemoteException
	{

		// Update the params based on new data
		double[][] result = update(currentDate);
		evaluate(currentDate);
	}

	/**
	 * To inform this instance to proceed to a new simulation day
	 * @param p_date The date of the new day
	 * @throws RemoteException This exception *MUST* be thrown by this method
	 */
	@Override
	public void proceedNewDay(int p_date)
		throws RemoteException
	{
		
		currentDate = p_date;
		double[][] params;

		// Calculate the parametes using the loop method
		params = calculateParams(p_date);
		double newPrice = calculatePrice(params);

		// If first day of the simulation run recursive method
		// to calculate all the old parameters
		if (currentDate == 101)
		{
			oldParams.clear();
			oldP.clear();
			oldParams.add(null);
			oldP.add(null);
			recursive(p_date-1);
		}

		// Else get the udated parameter
		else
		{
			// If simulation lasts more than one day update for the previous days data
			if (steps > 1)
			{
				double[][] result = update(currentDate-1);
			}
			params = oldParams.get(p_date-1);

		}

		double l_newPrice = calculatePrice(params);
		m_platformStub.publishPrice(m_type, (float)l_newPrice);

	}


	public static void main(final String[] p_args)
		throws RemoteException, NotBoundException
	{

		new PseudoLeader();
	}

	// Evaluate the quality of the solution
	// Calculates Root Mean Square Errors (rmse),
	// Mean Absolute Percentage Error (mape), and
	// R-square to evaluate the results
	public void evaluate(int day)
		throws RemoteException
	{
		Record record;
		double leaderPrice;
		double followerPrice;

		double expFollowerPrice;

		double[][] params;

		double rmse=0;
		double mape=0;
		double sse=0;
		double sst=0;
		double followePriceSum = 0;

		for(int i=1; i<=day; i++)
		{
			record = m_platformStub.query(PlayerType.LEADER, i);
			leaderPrice = record.m_leaderPrice;
			followerPrice = record.m_followerPrice;
			followePriceSum += followerPrice;

			params = oldParams.get(i);

			expFollowerPrice = params[0][0]+params[1][0]*leaderPrice;

			rmse += Math.pow(followerPrice - expFollowerPrice, 2);
			mape += Math.abs((followerPrice - expFollowerPrice) / followerPrice);
			sse += Math.pow(followerPrice - expFollowerPrice, 2);
		}


		double y = (1.0/(double)day) * followePriceSum;

		for(int i=1; i<=day; i++)
		{
			record = m_platformStub.query(PlayerType.LEADER, i);
			sst += Math.pow(record.m_followerPrice - y, 2);
		}

		System.out.println("RMSE: " + (Math.sqrt((1.0/(double)day)*rmse)));
		System.out.println("MAPE: " + ((1.0/(double)day)*mape));
		System.out.println("R-square: " + (1-(sse/sst)));
	}

	// Calculates new price based on the parameters
	public double calculatePrice(double[][] params)
	{
		return ((1-0.3*params[1][0])*COST+2+0.3*params[0][0])/(2-0.6*params[1][0]);
	}

	// Updates the parameters when new data is available
	// Returns the updated parameters
	public double[][] update(int day)
		throws RemoteException
	{
		Record record = m_platformStub.query(PlayerType.LEADER, day);
		double [][] leaderPrice = {{1.0, record.m_leaderPrice}};
		double predError = record.m_followerPrice - matrixMultiply(leaderPrice, oldParams.get(day-1))[0][0]; 
		double[][] adjFac = calculateAdjustingFactor(day, oldP.get(day-1));
		double[][] result = matrixAdd(oldParams.get(day-1), matrixScalarMultiply(adjFac, predError));
		
		if(oldParams.size() > day){
			oldParams.subList(day, oldParams.size()).clear();
		}
		oldParams.add(result);
		return result;

	}

	public double[][] calculateParams(int date)
		throws RemoteException
	{
		Record record;
		double leaderPrice;
		double followerPrice;
		double uLSquared = 0;
		double uL = 0;
		double uF = 0;
		double uLuF = 0;

		for(int i=1; i<=date; i++)
		{
			record = m_platformStub.query(PlayerType.LEADER, i);
			leaderPrice = record.m_leaderPrice;
			followerPrice = record.m_followerPrice;

			uLSquared += Math.pow(leaderPrice, 2);
			uL += leaderPrice;
			uF += followerPrice;
			uLuF += (leaderPrice*followerPrice);
		}


		double a = (uLSquared*uF - uL*uLuF) / (date*uLSquared - Math.pow(uL, 2));

		double b = (date*uLuF - uL*uF) / (date*uLSquared - Math.pow(uL, 2));

		double[][] result = {{a}, {b}};

		return result;
	}


	// Recursuve Least Square Approach
	// Calculates the parameters for all dates till parameter day
	// Parameters for each day stored in oldParams
	// Returns the most recent parameters
	public double[][] recursive(int day)
		throws RemoteException
	{
		if (day == 1)
		{
			Record record = m_platformStub.query(PlayerType.LEADER, day);
			double [][] leaderPrice = {{1.0}, {record.m_leaderPrice}};
			Pair initialPTheta = initialize(100);
			double predError = record.m_followerPrice - matrixMultiply(matrixTranspose(leaderPrice), initialPTheta.getThetaMatrix())[0][0]; 
			double[][] adjFac = calculateAdjustingFactor(day, initialPTheta.getPMatrix());
			double[][] result = matrixAdd(initialPTheta.getThetaMatrix(), matrixScalarMultiply(adjFac, predError));
			oldParams.add(result);
			return result;
		}

		double[][] oldParam = recursive(day-1);
		double[][] p = calculateP(day);
		double[][] adjFac = calculateAdjustingFactor(day, p);
		Record record = m_platformStub.query(PlayerType.LEADER, day);
		double[][] leaderPrice = {{1.0}, {record.m_leaderPrice}};
		double predError = record.m_followerPrice - matrixMultiply(matrixTranspose(leaderPrice), oldParam)[0][0]; 
		double[][] result = matrixAdd(oldParam, matrixScalarMultiply(adjFac, predError));
		oldParams.add(result);

		return result;
	}


	// Calculates the Adjusting Factor in RLS Approach
	public double[][] calculateAdjustingFactor(int day, double[][] p)
		throws RemoteException
	{
		Record record = m_platformStub.query(PlayerType.LEADER, day);
		double [][] leaderPrice = {{1.0}, {record.m_leaderPrice}};
		double [][] top = matrixMultiply(p, leaderPrice);
		double bottom = LAMBDA + matrixMultiply(matrixTranspose(leaderPrice), top)[0][0];
		return matrixScalarDivide(top, bottom);
	}

	// Recursively calculates the parameter P in RLS Approach
	// Stores the value of P for each day in oldP
	public double[][] calculateP(int day)
		throws RemoteException
	{
		if (day == 1)
		{
			if (oldP.get(0) == null)
			{
				return initialize(100).getPMatrix();
			}

			else
				return oldP.get(0);
		}	

		double[][] p = calculateP(day-1);
		Record record = m_platformStub.query(PlayerType.LEADER, day);
		double[][] leaderPrice = {{1.0}, {record.m_leaderPrice}};
		double[][] top = matrixMultiply(p, leaderPrice);
		double[][] top2 = matrixMultiply(matrixTranspose(leaderPrice), p);
		double[][] top3 = matrixMultiply(top, top2);
		double [][] bottom1 = matrixMultiply(p, leaderPrice);
		double bottom = LAMBDA + matrixMultiply(matrixTranspose(leaderPrice), bottom1)[0][0];
		double[][] divison = matrixScalarDivide(top3, bottom);
		double[][] result = matrixScalarMultiply(matrixSubtract(p, divison), (1/LAMBDA));
		oldP.add(result);
		return result;


	}

	// Calculates the initial values for P and theta (params)
	// Returns Pair containing the initial matrices
	public Pair initialize(int t)
		throws RemoteException
	{

		Record record;
		double[][] sumP = {{0.0f, 0.0f}, {0.0f, 0.0f}};
		double[][] sumTheta = {{0.0f}, {0.0f}};
		double[][] leaderPrice = new double[2][1];
		double[][] transpose = new double[1][2];
		double followerPrice;


		for(int i=1; i<=t; i++)
		{
			record = m_platformStub.query(PlayerType.LEADER, i);
			leaderPrice[0][0] = 1.0f;
			leaderPrice[1][0] = record.m_leaderPrice;
			transpose[0][0] = leaderPrice[0][0];
			transpose[0][1] = leaderPrice[1][0];
			followerPrice = record.m_followerPrice;
			sumP = matrixAdd(sumP, matrixScalarMultiply(matrixMultiply(leaderPrice, transpose), Math.pow(LAMBDA, (t-i))));
			sumTheta = matrixAdd(sumTheta, matrixScalarMultiply(matrixScalarMultiply(leaderPrice, followerPrice), Math.pow(LAMBDA, (t-i))));
		}

		double[][] finalTheta = matrixMultiply(matrixInverse(sumP), sumTheta);
		oldParams.set(0, finalTheta);
		oldP.set(0, sumP);

		return new Pair(sumP, finalTheta);

	}

	public static double[][] matrixTranspose(double [][] m){
        double[][] temp = new double[m[0].length][m.length];
        for (int i = 0; i < m.length; i++)
            for (int j = 0; j < m[0].length; j++)
                temp[j][i] = m[i][j];
        return temp;
    }


    public static double[][] matrixMultiply(double[][] m1, double[][] m2) {
        int m1ColLength = m1[0].length; // m1 columns length
        int m2RowLength = m2.length;    // m2 rows length
        if(m1ColLength != m2RowLength)
        {
        	System.out.println("matrix multiplication not possible");
        	return null; // matrix multiplication is not possible	
        } 
        int mRRowLength = m1.length;    // m result rows length
        int mRColLength = m2[0].length; // m result columns length
        double[][] mResult = new double[mRRowLength][mRColLength];
        for(int i = 0; i < mRRowLength; i++) {         // rows from m1
            for(int j = 0; j < mRColLength; j++) {     // columns from m2
                for(int k = 0; k < m1ColLength; k++) { // columns from m1
                    mResult[i][j] += m1[i][k] * m2[k][j];
                }
            }
        }
        return mResult;
    }


    public static double[][] matrixScalarMultiply(double[][] m, double f)
    {
    	double[][] temp = new double[m.length][m[0].length];
    	for (int i = 0; i < m.length; i++)
            for (int j = 0; j < m[0].length; j++)
                temp[i][j] = m[i][j] * f;
        return temp;
    }

    public static double[][] matrixScalarDivide(double[][] m, double f)
    {
    	double[][] temp = new double[m.length][m[0].length];
    	for (int i = 0; i < m.length; i++)
            for (int j = 0; j < m[0].length; j++)
                temp[i][j] = m[i][j] / f;
        return temp;
    }

    public static double[][] matrixAdd(double[][] m1, double[][] m2)
    {
    	double[][] temp = new double[m1.length][m1[0].length];
    	for (int i = 0; i < m1.length; i++)
            for (int j = 0; j < m1[0].length; j++)
                temp[i][j] = m1[i][j] + m2[i][j];
        return temp;
    }

    public static double[][] matrixSubtract(double[][] m1, double[][] m2)
    {
    	double[][] temp = new double[m1.length][m1[0].length];
    	for (int i = 0; i < m1.length; i++)
            for (int j = 0; j < m1[0].length; j++)
                temp[i][j] = m1[i][j] - m2[i][j];
        return temp;
    }

    public static double[][] matrixInverse(double[][] m)
    {
    	double[][] temp = new double[m.length][m[0].length];
    	temp[0][0] = m[1][1];
    	temp[0][1] = -m[0][1];
    	temp[1][0] = -m[1][0];
    	temp[1][1] = m[0][0];
    	
    	double s = 1.0 / (m[0][0]*m[1][1] - m[0][1]*m[1][0]);
    	
    	temp = matrixScalarMultiply(temp, s);

    	return temp;
    }

    public static void printMatrix(double[][] matrix) {
        for (double[] array : matrix) {
        	for (double i : array) {
	            System.out.print(i);
    	        System.out.print("\t");
        	}
        	System.out.println();
        }
        
    }

    public class Pair
    {
    	private double[][] p;
    	private double[][] theta;

    	public Pair(double[][] p, double[][] theta)
    	{
    		this.p = p;
    		this.theta = theta;
    	}

    	public double[][] getPMatrix()
    	{
    		return p;
    	}

    	public double[][] getThetaMatrix()
    	{
    		return theta;
    	}
    }

}
