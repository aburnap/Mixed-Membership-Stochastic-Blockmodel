// (C) Copyright 2009 Microsoft Research Cambridge
using System;
using System.Collections.Generic;
using System.Text;
using MicrosoftResearch.Infer.Distributions;
using MicrosoftResearch.Infer.Maths;

namespace MicrosoftResearch.Infer.Factors
{
	/// <summary>
	/// Provides outgoing messages for <see cref="Factor.DiscreteUniform"/>, given random arguments to the function.
	/// </summary>
	/// <remarks>Factor is f(sample, size) = 1(sample &lt; size)/size</remarks>
	[FactorMethod(typeof(Factor), "DiscreteUniform")]
	[Quality(QualityBand.Stable)]
	public static class DiscreteUniform
	{
		/// <summary>
		/// Evidence message for EP
		/// </summary>
		/// <param name="sample">Constant value for 'sample'.</param>
		/// <param name="size">Constant value for 'size'.</param>
		/// <returns>Logarithm of the factor's average value across the given argument distributions</returns>
		/// <remarks><para>
		/// The formula for the result is <c>log(factor(sample,size))</c>.
		/// </para></remarks>
		public static double LogAverageFactor(int sample, int size)
		{
			return (sample < size) ? -Math.Log(size) : Double.NegativeInfinity;
		}
		/// <summary>
		/// Evidence message for EP
		/// </summary>
		/// <param name="sample">Constant value for 'sample'.</param>
		/// <param name="size">Constant value for 'size'.</param>
		/// <returns>Logarithm of the factor's contribution the EP model evidence</returns>
		/// <remarks><para>
		/// The formula for the result is <c>log(factor(sample,size))</c>.
		/// Adding up these values across all factors and variables gives the log-evidence estimate for EP.
		/// </para></remarks>
		public static double LogEvidenceRatio(int sample, int size) { return LogAverageFactor(sample, size); }

		/// <summary>
		/// Evidence message for EP
		/// </summary>
		/// <param name="sample">Constant value for 'sample'.</param>
		/// <param name="size">Incoming message from 'size'.</param>
		/// <returns>Logarithm of the factor's average value across the given argument distributions</returns>
		/// <remarks><para>
		/// The formula for the result is <c>log(sum_(size) p(size) factor(sample,size))</c>.
		/// </para></remarks>
		public static double LogAverageFactor(int sample, Discrete size)
		{
			double z = 0.0;
			for (int i = sample+1; i < size.Dimension; i++)
			{
				z += size[i]/i;
			}
			return Math.Log(z);
		}
		/// <summary>
		/// Evidence message for EP
		/// </summary>
		/// <param name="sample">Constant value for 'sample'.</param>
		/// <param name="size">Incoming message from 'size'.</param>
		/// <returns>Logarithm of the factor's contribution the EP model evidence</returns>
		/// <remarks><para>
		/// The formula for the result is <c>log(sum_(size) p(size) factor(sample,size))</c>.
		/// Adding up these values across all factors and variables gives the log-evidence estimate for EP.
		/// </para></remarks>
		public static double LogEvidenceRatio(int sample, Discrete size)
		{
			return LogAverageFactor(sample, size);
		}

		/// <summary>
		/// Evidence message for EP
		/// </summary>
		/// <param name="sample">Incoming message from 'sample'.</param>
		/// <param name="to_sample">Outgoing message to 'sample'.</param>
		/// <returns>Logarithm of the factor's average value across the given argument distributions</returns>
		/// <remarks><para>
		/// The formula for the result is <c>log(sum_(sample) p(sample) factor(sample,size))</c>.
		/// </para></remarks>
		public static double LogAverageFactor(Discrete sample, [Fresh] Discrete to_sample)
		{
			return to_sample.GetLogAverageOf(sample);
		}
		/// <summary>
		/// Evidence message for EP
		/// </summary>
		/// <param name="sample">Incoming message from 'sample'. Must be a proper distribution.  If uniform, the result will be uniform.</param>
		/// <param name="size">Incoming message from 'size'.</param>
		/// <returns>Logarithm of the factor's contribution the EP model evidence</returns>
		/// <remarks><para>
		/// The formula for the result is <c>log(sum_(sample,size) p(sample,size) factor(sample,size) / sum_sample p(sample) messageTo(sample))</c>.
		/// Adding up these values across all factors and variables gives the log-evidence estimate for EP.
		/// </para></remarks>
		/// <exception cref="ImproperMessageException"><paramref name="sample"/> is not a proper distribution</exception>
		public static double LogEvidenceRatio([SkipIfUniform] Discrete sample, Discrete size)
		{
			return Math.Log(1-size[0]);
		}
		/// <summary>
		/// Evidence message for EP
		/// </summary>
		/// <param name="sample">Incoming message from 'sample'. Must be a proper distribution.  If uniform, the result will be uniform.</param>
		/// <param name="size">Constant value for 'size'.</param>
		/// <returns>Logarithm of the factor's contribution the EP model evidence</returns>
		/// <remarks><para>
		/// The formula for the result is <c>log(sum_(sample) p(sample) factor(sample,size) / sum_sample p(sample) messageTo(sample))</c>.
		/// Adding up these values across all factors and variables gives the log-evidence estimate for EP.
		/// </para></remarks>
		/// <exception cref="ImproperMessageException"><paramref name="sample"/> is not a proper distribution</exception>
		public static double LogEvidenceRatio([SkipIfUniform] Discrete sample, int size)
		{
			if (size == 0)
				return Double.NegativeInfinity;
			else
				return 0.0;
		}

		/// <summary>
		/// EP message to 'sample'
		/// </summary>
		/// <param name="size">Constant value for 'size'.</param>
		/// <param name="result">Modified to contain the outgoing message</param>
		/// <returns><paramref name="result"/></returns>
		/// <remarks><para>
		/// The outgoing message is the factor viewed as a function of 'sample' conditioned on the given values.
		/// </para></remarks>
		public static Discrete SampleAverageConditional(int size, Discrete result)
		{
			if (size == 0) return result; //throw new AllZeroException();
			if (result.Dimension < size) throw new ArgumentException("result.Dimension ("+result.Dimension+") < size ("+size+")");
			Vector probs = result.GetWorkspace();
			double invSize = 1.0/size;
			for (int i = 0; i < size; i++)
			{
				probs[i] = invSize;
			}
			for (int i = size; i < probs.Count; i++)
			{
				probs[i] = 0.0;
			}
			result.SetProbs(probs);
			return result;
		}

		/// <summary>
		/// EP message to 'sample'
		/// </summary>
		/// <param name="size">Incoming message from 'size'.</param>
		/// <param name="result">Modified to contain the outgoing message</param>
		/// <returns><paramref name="result"/></returns>
		/// <remarks><para>
		/// The outgoing message is a distribution matching the moments of 'sample' as the random arguments are varied.
		/// The formula is <c>proj[p(sample) sum_(size) p(size) factor(sample,size)]/p(sample)</c>.
		/// </para></remarks>
		public static Discrete SampleAverageConditional(Discrete size, Discrete result)
		{
			if (size.IsPointMass) return SampleAverageConditional(size.Point, result);
			if (result.Dimension < size.Dimension-1) throw new ArgumentException("result.Dimension ("+result.Dimension+") < size.Dimension-1 ("+size.Dimension+"-1)");
			Vector probs = result.GetWorkspace();
			// if size.Dimension = 8, then size ranges 0..7, which means sample ranges 0..6
			for (int i = probs.Count - 1; i >= size.Dimension-1; i--)
			{
				probs[i] = 0.0;
			}
			if (size.Dimension > 1)
			{
				// p(sample) = sum_size 1(sample < size)/size p(size)
				//           = sum_(size > sample) p(size)/size
				// e.g. if dimension == 4 then p(sample=0) = p(size=1)/1 + p(size=2)/2 + p(size=3)/3
				// note the value of size[0] is irrelevant
				probs[size.Dimension-2] = size[size.Dimension-1]/(size.Dimension-1);
				for (int i = size.Dimension - 3; i >= 0; i--)
				{
					probs[i] = probs[i+1] + size[i+1]/(i+1);
				}
			}
			result.SetProbs(probs);
			return result;
		}

		[Skip]
		public static Discrete SampleAverageConditionalInit(int size)
		{
			return Discrete.Uniform(size);
		}
		[Skip]
		public static Discrete SampleAverageConditionalInit([IgnoreDependency] Discrete size)
		{
			// if size.Dimension = 8, then size ranges 0..7, which means sample ranges 0..6
			return Discrete.Uniform(size.Dimension-1);
		}

		/// <summary>
		/// EP message to 'size'
		/// </summary>
		/// <param name="sample">Constant value for 'sample'.</param>
		/// <param name="result">Modified to contain the outgoing message</param>
		/// <returns><paramref name="result"/></returns>
		/// <remarks><para>
		/// The outgoing message is the factor viewed as a function of 'size' conditioned on the given values.
		/// </para></remarks>
		public static Discrete SizeAverageConditional(int sample, Discrete result)
		{
			Vector probs = result.GetWorkspace();
			for (int size = 0; size <= sample; size++)
			{
				probs[size] = 0.0;
			}
			for (int size = sample+1; size < probs.Count; size++)
			{
				probs[size] = 1.0/size;
			}
			result.SetProbs(probs);
			return result;
		}

		/// <summary>
		/// EP message to 'size'
		/// </summary>
		/// <param name="sample">Incoming message from 'sample'. Must be a proper distribution.  If uniform, the result will be uniform.</param>
		/// <param name="result">Modified to contain the outgoing message</param>
		/// <returns><paramref name="result"/></returns>
		/// <remarks><para>
		/// The outgoing message is a distribution matching the moments of 'size' as the random arguments are varied.
		/// The formula is <c>proj[p(size) sum_(sample) p(sample) factor(sample,size)]/p(size)</c>.
		/// </para></remarks>
		/// <exception cref="ImproperMessageException"><paramref name="sample"/> is not a proper distribution</exception>
		public static Discrete SizeAverageConditional([SkipIfUniform] Discrete sample, Discrete result)
		{
			if (sample.IsPointMass) return SizeAverageConditional(sample.Point, result);
			// p(size) =propto sum_sample p(sample) 1(sample < size)/size
			//         =propto p(sample < size)/size
			// e.g. if dimension == 4 then p(size=3) = p(sample < 3)/3
			Vector probs = result.GetWorkspace();
			probs[0] = 0.0;
			for (int size = 1; size < probs.Count; size++)
			{
				probs[size] = probs[size-1] + sample[size-1];
			}
			for (int size = 1; size < probs.Count; size++)
			{
				probs[size] /= size;
			}
			result.SetProbs(probs);
			return result;
		}

		// VMP -----------------------------------------------------------------

		/// <summary>
		/// Evidence message for VMP
		/// </summary>
		/// <param name="sample">Incoming message from 'sample'.</param>
		/// <param name="size">Constant value for 'size'.</param>
		/// <returns>Average of the factor's log-value across the given argument distributions</returns>
		/// <remarks><para>
		/// The formula for the result is <c>sum_(sample) p(sample) log(factor(sample,size))</c>.
		/// Adding up these values across all factors and variables gives the log-evidence estimate for VMP.
		/// </para></remarks>
		public static double AverageLogFactor(Discrete sample, int size)
		{
			double sum = 0.0;
			for (int i = 0; (i < sample.Dimension) && (i < size); i++)
			{
				sum += sample[i];
			}
			return -Math.Log(size)*sum;
		}
		/// <summary>
		/// Evidence message for VMP
		/// </summary>
		/// <param name="sample">Constant value for 'sample'.</param>
		/// <param name="size">Incoming message from 'size'.</param>
		/// <returns>Average of the factor's log-value across the given argument distributions</returns>
		/// <remarks><para>
		/// The formula for the result is <c>sum_(size) p(size) log(factor(sample,size))</c>.
		/// Adding up these values across all factors and variables gives the log-evidence estimate for VMP.
		/// </para></remarks>
		public static double AverageLogFactor(int sample, Discrete size)
		{
			double logZ = Double.NegativeInfinity;
			for (int n = 1; n < size.Dimension; n++) {
				logZ = MMath.LogSumExp(logZ, AverageLogFactor(sample, n));
			}
			return logZ;
		}
		/// <summary>
		/// Evidence message for VMP
		/// </summary>
		/// <param name="sample">Constant value for 'sample'.</param>
		/// <param name="size">Constant value for 'size'.</param>
		/// <returns>Average of the factor's log-value across the given argument distributions</returns>
		/// <remarks><para>
		/// The formula for the result is <c>log(factor(sample,size))</c>.
		/// Adding up these values across all factors and variables gives the log-evidence estimate for VMP.
		/// </para></remarks>
		public static double AverageLogFactor(int sample, int size)
		{
			return LogAverageFactor(sample, size);
		}
		/// <summary>
		/// Evidence message for VMP
		/// </summary>
		/// <param name="sample">Incoming message from 'sample'.</param>
		/// <param name="size">Incoming message from 'size'.</param>
		/// <returns>Average of the factor's log-value across the given argument distributions</returns>
		/// <remarks><para>
		/// The formula for the result is <c>sum_(sample,size) p(sample,size) log(factor(sample,size))</c>.
		/// Adding up these values across all factors and variables gives the log-evidence estimate for VMP.
		/// </para></remarks>
		public static double AverageLogFactor(Discrete sample, Discrete size)
		{
			double logZ = Double.NegativeInfinity;
			for (int n = 1; n < size.Dimension; n++)
			{
				logZ = MMath.LogSumExp(logZ, AverageLogFactor(sample, n));
			}
			return logZ;
		}

		[Skip]
		public static Discrete SampleAverageLogarithmInit(int size)
		{
			return Discrete.Uniform(size);
		}
		[Skip]
		public static Discrete SampleAverageLogarithmInit([IgnoreDependency] Discrete size)
		{
			// if size.Dimension = 8, then size ranges 0..7, which means sample ranges 0..6
			return Discrete.Uniform(size.Dimension-1);
		}
	
		/// <summary>
		/// VMP message to 'sample'
		/// </summary>
		/// <param name="size">Constant value for 'size'.</param>
		/// <param name="result">Modified to contain the outgoing message</param>
		/// <returns><paramref name="result"/></returns>
		/// <remarks><para>
		/// The outgoing message is the factor viewed as a function of 'sample' conditioned on the given values.
		/// </para></remarks>
		public static Discrete SampleAverageLogarithm(int size, Discrete result)
		{
			return SampleAverageConditional(size, result);
		}
		/// <summary>
		/// VMP message to 'sample'
		/// </summary>
		/// <param name="size">Incoming message from 'size'.</param>
		/// <param name="result">Modified to contain the outgoing message</param>
		/// <returns><paramref name="result"/></returns>
		/// <remarks><para>
		/// The outgoing message is the exponential of the average log-factor value, where the average is over all arguments except 'sample'.
		/// The formula is <c>exp(sum_(size) p(size) log(factor(sample,size)))</c>.
		/// </para></remarks>
		public static Discrete SampleAverageLogarithm(Discrete size, Discrete result)
		{
			return SampleAverageConditional(size, result);
		}
		/// <summary>
		/// VMP message to 'size'
		/// </summary>
		/// <param name="sample">Constant value for 'sample'.</param>
		/// <param name="result">Modified to contain the outgoing message</param>
		/// <returns><paramref name="result"/></returns>
		/// <remarks><para>
		/// The outgoing message is the factor viewed as a function of 'size' conditioned on the given values.
		/// </para></remarks>
		public static Discrete SizeAverageLogarithm(int sample, Discrete result)
		{
			return SizeAverageConditional(sample, result);
		}
		/// <summary>
		/// VMP message to 'size'
		/// </summary>
		/// <param name="sample">Incoming message from 'sample'. Must be a proper distribution.  If uniform, the result will be uniform.</param>
		/// <param name="result">Modified to contain the outgoing message</param>
		/// <returns><paramref name="result"/></returns>
		/// <remarks><para>
		/// The outgoing message is the exponential of the average log-factor value, where the average is over all arguments except 'size'.
		/// The formula is <c>exp(sum_(sample) p(sample) log(factor(sample,size)))</c>.
		/// </para></remarks>
		/// <exception cref="ImproperMessageException"><paramref name="sample"/> is not a proper distribution</exception>
		public static Discrete SizeAverageLogarithm([SkipIfUniform] Discrete sample, Discrete result)
		{
			return SizeAverageConditional(sample, result);
		}
	}
}
