// (C) Copyright 2008 Microsoft Research Cambridge
using System;
using System.Collections.Generic;
using System.Text;
using MicrosoftResearch.Infer.Maths;
using MicrosoftResearch.Infer.Collections;
using MicrosoftResearch.Infer.Utils;
using MicrosoftResearch.Infer.Factors;

namespace MicrosoftResearch.Infer.Distributions
{
	/// <summary>
	/// Binomial distribution over the integers [0,n]
	/// </summary>
	/// <remarks>
	/// The formula for the distribution is <c>p(x) = n!/(n-x)!/x! p^x (1-p)^(n-x)</c>.
	/// In this implementation, we use a generalization that includes two extra shape parameters (a,b)
	/// The formula for the generalized distribution is <c>p(x) =propto 1/x!^a 1/(n-x)!^b exp(x*logOdds)</c>.
	/// With this extension, we can represent a uniform distribution via (logOdds=0,a=0,b=0) and
	/// a point mass via logOdds=+/-infinity or a=infinity or b=infinity.
	/// This family is closed under multiplication, while the standard Binomial is not.
	/// </remarks>
	[Serializable]
	[Quality(QualityBand.Experimental)]
	public struct Binomial : IDistribution<int>,
		SettableTo<Binomial>, SettableToProduct<Binomial>, SettableToRatio<Binomial>,
		SettableToPower<Binomial>, SettableToWeightedSum<Binomial>,
		CanGetLogAverageOf<Binomial>,
		CanGetAverageLog<Binomial>,
		Sampleable<int>, CanGetMean<double>, CanGetVariance<double>,
		CanGetMeanAndVarianceOut<double, double>
	{
		public int TrialCount;
		public double LogOdds;
		public double A, B;

		public double ProbSuccess { get { return MMath.Logistic(LogOdds); } }
		public double ProbFailure { get { return MMath.Logistic(-LogOdds); } }

		public Binomial(int trialCount, double probSuccess)
		{
			TrialCount = trialCount;
			LogOdds = MMath.Logit(probSuccess);
			A = 1;
			B = 1;
		}

		[Construction("TrialCount", "LogOdds", "A", "B")]
		public static Binomial FromNatural(int trialCount, double logOdds, double a = 1, double b = 1)
		{
			Binomial result = new Binomial();
			result.TrialCount = trialCount;
			result.LogOdds = logOdds;
			result.A = a;
			result.B = b;
			return result;
		}

		public static Binomial Uniform(int trialCount)
		{
			return Binomial.FromNatural(trialCount, 0, 0, 0);
		}

		public void SetTo(Binomial that)
		{
			this.TrialCount = that.TrialCount;
			this.LogOdds = that.LogOdds;
			this.A = that.A;
			this.B = that.B;
		}

		public double GetMean()
		{
			if (A != 1 || B != 1) throw new NotImplementedException("A != 1 or B !=1 not implemented");
			return TrialCount*ProbSuccess;
		}

		public double GetVariance()
		{
			if (A != 1 || B != 1) throw new NotImplementedException("A != 1 or B !=1 not implemented");
			return TrialCount*ProbSuccess*ProbFailure;
		}

		public void GetMeanAndVariance(out double mean, out double variance)
		{
			if (A != 1 || B != 1) throw new NotImplementedException("A != 1 or B !=1 not implemented");
			mean = TrialCount*ProbSuccess;
			variance = mean*ProbFailure;
		}

		public void SetToProduct(Binomial a, Binomial b)
		{
			if (a.TrialCount != b.TrialCount) throw new ArgumentException("a.TrialCount ("+a.TrialCount+") != b.TrialCount ("+b.TrialCount+")");
			throw new NotImplementedException();
			TrialCount = a.TrialCount;
			LogOdds = a.LogOdds + b.LogOdds;
			A = a.A + b.A;
			B = a.B + b.B;
		}

		public void SetToRatio(Binomial numerator, Binomial denominator)
		{
			if (numerator.TrialCount != denominator.TrialCount) throw new ArgumentException("numerator.TrialCount ("+numerator.TrialCount+") != denominator.TrialCount ("+denominator.TrialCount+")");
			throw new NotImplementedException();
			TrialCount = numerator.TrialCount;
			LogOdds = numerator.LogOdds - denominator.LogOdds;
			A = numerator.A - denominator.A;
			B = numerator.B - denominator.B;
		}

		public void SetToPower(Binomial value, double exponent)
		{
			throw new NotImplementedException();
			LogOdds = exponent*LogOdds;
			A = exponent*A;
			B = exponent*B;
		}

		public void SetToSum(double weight1, Binomial value1, double weight2, Binomial value2)
		{
			throw new NotImplementedException();
		}

		public double GetLogAverageOf(Binomial that)
		{
			throw new NotImplementedException();
		}

		public double GetAverageLog(Binomial that)
		{
			throw new NotImplementedException();
		}

		public int Sample()
		{
			if(A != 1 || B != 1) throw new NotImplementedException();
			return Rand.Binomial(TrialCount, ProbSuccess);
		}

		public int Sample(int result)
		{
			return Sample();
		}

		public object Clone()
		{
			Binomial result = new Binomial();
			result.SetTo(this);
			return result;
		}

		public int Point
		{
			get
			{
				throw new NotImplementedException();
			}
			set
			{
				throw new NotImplementedException();
			}
		}

		public bool IsPointMass
		{
			get { return double.IsInfinity(LogOdds) || double.IsInfinity(A) || double.IsInfinity(B); }
		}

		public double MaxDiff(object that)
		{
			if(!(that is Binomial)) return double.PositiveInfinity;
			Binomial thatd = (Binomial)that;
			return Math.Max(Math.Abs(TrialCount - thatd.TrialCount), Math.Abs(LogOdds - thatd.LogOdds));
		}

		public void SetToUniform()
		{
			LogOdds = 0;
			A = 0;
			B = 0;
		}

		public bool IsUniform()
		{
			return (LogOdds == 0) && (A == 0) && (B == 0);
		}

		public double GetLogProb(int value)
		{
			if (value < 0 || value > TrialCount) return double.NegativeInfinity;
			return MMath.GammaLn(TrialCount+1) - A*MMath.GammaLn(value+1) - B*MMath.GammaLn(TrialCount-value+1) + value*LogOdds + TrialCount*MMath.LogisticLn(-LogOdds);
		}
	}
}
