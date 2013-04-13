// (C) Copyright 2008 Microsoft Research Cambridge
using System;
using System.Collections.Generic;
using System.Text;
using MicrosoftResearch.Infer.Distributions;
using MicrosoftResearch.Infer;
using MicrosoftResearch.Infer.Maths;

namespace MicrosoftResearch.Infer.Factors
{
	/// <summary>
	/// Provides outgoing messages for <see cref="Factor.Copy{T}"/>, given random arguments to the function.
	/// </summary>
	[FactorMethod(typeof(Factor), "Copy<>")]
	[Quality(QualityBand.Mature)]
	public static class CopyOp<T>
	{
		/// <summary>
		/// Evidence message for EP
		/// </summary>
		/// <returns>Logarithm of the factor's average value across the given argument distributions</returns>
		/// <remarks><para>
		/// The formula for the result is <c>log(factor(copy,value))</c>.
		/// </para></remarks>
		public static double LogAverageFactor(T copy, T value)
		{
			IEqualityComparer<T> equalityComparer = Utils.Util.GetEqualityComparer<T>();
			return equalityComparer.Equals(copy, value) ? 0.0 : Double.NegativeInfinity;
		}
		public static double LogAverageFactor<TDist>(TDist copy, T value)
			where TDist : CanGetLogProb<T>
		{
			return copy.GetLogProb(value);
		}
		public static double LogAverageFactor<TDist>(TDist copy, TDist value)
			where TDist : CanGetLogAverageOf<TDist>
		{
			return value.GetLogAverageOf(copy);
		}
		public static double LogAverageFactor<TDist>(T copy, TDist value)
			where TDist : CanGetLogProb<T>
		{
			return value.GetLogProb(copy);
		}

		public static double LogEvidenceRatio(T copy, T value) { return LogAverageFactor(copy, value); }
		/// <summary>
		/// Evidence message for EP
		/// </summary>
		/// <param name="copy">Incoming message from 'copy'.</param>
		/// <returns>Logarithm of the factor's contribution the EP model evidence</returns>
		/// <remarks><para>
		/// The formula for the result is <c>log(sum_(copy) p(copy) factor(copy,value) / sum_copy p(copy) messageTo(copy))</c>.
		/// Adding up these values across all factors and variables gives the log-evidence estimate for EP.
		/// </para></remarks>
		[Skip]
		public static double LogEvidenceRatio<TDist>(TDist copy)
			where TDist : IDistribution<T>
		{
			return 0.0;
		}

		/// <summary>
		/// Evidence message for EP
		/// </summary>
		/// <param name="copy">Constant value for 'copy'.</param>
		/// <param name="value">Incoming message from 'value'.</param>
		/// <returns>Logarithm of the factor's contribution the EP model evidence</returns>
		/// <remarks><para>
		/// The formula for the result is <c>log(sum_(value) p(value) factor(copy,value))</c>.
		/// Adding up these values across all factors and variables gives the log-evidence estimate for EP.
		/// </para></remarks>
		public static double LogEvidenceRatio<TDist>(T copy, TDist value)
			where TDist : CanGetLogProb<T>
		{
			return value.GetLogProb(copy);
		}

		/// <summary>
		/// EP message to 'value'
		/// </summary>
		/// <param name="copy">Incoming message from 'copy'.</param>
		/// <returns>The outgoing EP message to the 'value' argument</returns>
		/// <remarks><para>
		/// The outgoing message is a distribution matching the moments of 'value' as the random arguments are varied.
		/// The formula is <c>proj[p(value) sum_(copy) p(copy) factor(copy,value)]/p(value)</c>.
		/// </para></remarks>
		public static TDist ValueAverageConditional<TDist>([IsReturned] TDist copy)
			where TDist : IDistribution<T>
		{
			return copy;
		}

		/// <summary>
		/// EP message to 'value'
		/// </summary>
		/// <param name="copy">Constant value for 'copy'.</param>
		/// <returns>The outgoing EP message to the 'value' argument</returns>
		/// <remarks><para>
		/// The outgoing message is the factor viewed as a function of 'value' conditioned on the given values.
		/// </para></remarks>
		public static T ValueAverageConditional([IsReturned] T copy)
		{
			return copy;
		}

		/// <summary>
		/// EP message to 'value'
		/// </summary>
		/// <param name="copy">Constant value for 'copy'.</param>
		/// <param name="result">Modified to contain the outgoing message</param>
		/// <returns><paramref name="result"/></returns>
		/// <remarks><para>
		/// The outgoing message is the factor viewed as a function of 'value' conditioned on the given values.
		/// </para></remarks>
		public static TDist ValueAverageConditional<TDist>(T copy, TDist result)
			where TDist : IDistribution<T>
		{
			result.Point = copy;
			return result;
		}

		/// <summary>
		/// EP message to 'copy'
		/// </summary>
		/// <param name="Value">Incoming message from 'value'.</param>
		/// <returns>The outgoing EP message to the 'copy' argument</returns>
		/// <remarks><para>
		/// The outgoing message is a distribution matching the moments of 'copy' as the random arguments are varied.
		/// The formula is <c>proj[p(copy) sum_(value) p(value) factor(copy,value)]/p(copy)</c>.
		/// </para></remarks>
		public static TDist CopyAverageConditional<TDist>([IsReturned] TDist Value)
			where TDist : IDistribution<T>
		{
			return Value;
		}

		/// <summary>
		/// EP message to 'copy'
		/// </summary>
		/// <param name="Value">Constant value for 'value'.</param>
		/// <returns>The outgoing EP message to the 'copy' argument</returns>
		/// <remarks><para>
		/// The outgoing message is the factor viewed as a function of 'copy' conditioned on the given values.
		/// </para></remarks>
		public static T CopyAverageConditional([IsReturned] T Value)
		{
			return Value;
		}

		//-- VMP ------------------------------------------------------------------------------------
		// NOTE: AverageLogFactor operators are explicit here so that the correct overload is used if
		// the copy or value is a truncated Gaussian

		/// <summary>
		/// Evidence message for VMP
		/// </summary>
		/// <returns>Zero</returns>
		/// <remarks><para>
		/// The formula for the result is <c>log(factor(copy,value))</c>.
		/// Adding up these values across all factors and variables gives the log-evidence estimate for VMP.
		/// </para></remarks>
		[Skip]
		public static double AverageLogFactor<TDist>(TDist copy, TDist Value)
				where TDist : IDistribution<T>
		{
			return 0.0;
		}

		/// <summary>
		/// Evidence message for VMP
		/// </summary>
		/// <returns>Zero</returns>
		/// <remarks><para>
		/// The formula for the result is <c>log(factor(copy,value))</c>.
		/// Adding up these values across all factors and variables gives the log-evidence estimate for VMP.
		/// </para></remarks>
		[Skip]
		public static double AverageLogFactor<TDist>(TDist copy, T Value)
				where TDist : IDistribution<T>
		{
			return 0.0;
		}

		/// <summary>
		/// Evidence message for VMP
		/// </summary>
		/// <returns>Zero</returns>
		/// <remarks><para>
		/// The formula for the result is <c>log(factor(copy,value))</c>.
		/// Adding up these values across all factors and variables gives the log-evidence estimate for VMP.
		/// </para></remarks>
		[Skip]
		public static double AverageLogFactor<TDist>(T copy, TDist Value)
				where TDist : IDistribution<T>
		{
			return 0.0;
		}

		public static double AverageLogFactor(T copy, T Value)
		{
			return LogAverageFactor(copy, Value);
		}

		// must have upward Trigger to match the Trigger on UsesEqualDef.UsesAverageLogarithm
		/// <summary>
		/// VMP message to 'value'
		/// </summary>
		/// <param name="copy">Incoming message from 'copy'.</param>
		/// <returns>The outgoing VMP message to the 'value' argument</returns>
		/// <remarks><para>
		/// The outgoing message is the factor viewed as a function of 'value' with 'copy' integrated out.
		/// The formula is <c>sum_copy p(copy) factor(copy,value)</c>.
		/// </para></remarks>
		public static TDist ValueAverageLogarithm<TDist>([IsReturned] TDist copy)
			where TDist : IDistribution<T>
		{
			return copy;
		}

		/// <summary>
		/// VMP message to 'value'
		/// </summary>
		/// <param name="copy">Constant value for 'copy'.</param>
		/// <param name="result">Modified to contain the outgoing message</param>
		/// <returns><paramref name="result"/></returns>
		/// <remarks><para>
		/// The outgoing message is the factor viewed as a function of 'value' conditioned on the given values.
		/// </para></remarks>
		public static TDist ValueAverageLogarithm<TDist>(T copy, TDist result)
			where TDist : IDistribution<T>
		{
			result.Point = copy;
			return result;
		}

		/// <summary>
		/// VMP message to 'copy'
		/// </summary>
		/// <param name="Value">Incoming message from 'value'.</param>
		/// <returns>The outgoing VMP message to the 'copy' argument</returns>
		/// <remarks><para>
		/// The outgoing message is a distribution matching the moments of 'copy' as the random arguments are varied.
		/// The formula is <c>proj[sum_(value) p(value) factor(copy,value)]</c>.
		/// </para></remarks>
		public static TDist CopyAverageLogarithm<TDist>([IsReturned] TDist Value)
			where TDist : IDistribution<T>
		{
			return Value;
		}

		[Skip]
		public static TDist CopyDeriv<TDist>(TDist result)
			where TDist : IDistribution<T>, SettableToUniform
		{
			result.SetToUniform();
			return result;
		}

		/* Methods for using the Copy factor to convert between Gaussians and Truncated Gaussians */

		/// <summary>
		/// VMP message to 'value'
		/// </summary>
		/// <param name="copy">Incoming message from 'copy'.</param>
		/// <returns>The outgoing VMP message to the 'value' argument</returns>
		/// <remarks><para>
		/// The outgoing message is the factor viewed as a function of 'value' with 'copy' integrated out.
		/// The formula is <c>sum_copy p(copy) factor(copy,value)</c>.
		/// </para></remarks>
		[Quality(QualityBand.Preview)]
		public static TruncatedGaussian ValueAverageLogarithm(Gaussian copy)
		{
			return new TruncatedGaussian(copy);
		}

		/// <summary>
		/// VMP message to 'copy'
		/// </summary>
		/// <param name="value">Incoming message from 'value'.</param>
		/// <returns>The outgoing VMP message to the 'copy' argument</returns>
		/// <remarks><para>
		/// The outgoing message is a distribution matching the moments of 'copy' as the random arguments are varied.
		/// The formula is <c>proj[sum_(value) p(value) factor(copy,value)]</c>.
		/// </para></remarks>
		[Quality(QualityBand.Preview)]
		public static Gaussian CopyAverageLogarithm(TruncatedGaussian value)
		{
			return value.ToGaussian();
		}

		/// <summary>
		/// VMP message to 'value'
		/// </summary>
		/// <param name="copy">Incoming message from 'copy'.</param>
		/// <param name="value">Incoming message from 'value'.</param>
		/// <param name="to_value">Previous outgoing message to 'value'.</param>
		/// <returns>The outgoing VMP message to the 'value' argument</returns>
		/// <remarks><para>
		/// The outgoing message is the factor viewed as a function of 'value' with 'copy' integrated out.
		/// The formula is <c>sum_copy p(copy) factor(copy,value)</c>.
		/// </para></remarks>
		[Quality(QualityBand.Preview)]
		public static Gaussian ValueAverageLogarithm(TruncatedGaussian copy, [Proper] Gaussian value, Gaussian to_value)
		{
			var a = value / to_value;
			copy *= new TruncatedGaussian(a); // is this ok? 
			var result = copy.ToGaussian() / a;
			return result;
		}

		/// <summary>
		/// Evidence message for VMP
		/// </summary>
		/// <param name="copy">Message from copy</param>
		/// <param name="value">Message from value</param>
		/// <param name="to_value">Message to value</param>
		/// In Variational Message Passing, the evidence contribution of a deterministic factor is zero.
		/// Adding up these values across all factors and variables gives the log-evidence estimate for VMP.
		/// This factor is implicitly maintaining the truncated Gaussian variational posterior. Therefore
		/// we need to remove the entropy of the Gaussian representation, and add the entropy for the
		/// truncated Gaussian
		/// <returns>Evidence</returns>
		[Quality(QualityBand.Preview)]
		public static double AverageLogFactor(TruncatedGaussian copy, Gaussian value, Gaussian to_value)
		{
			var a = value / to_value;
			copy *= new TruncatedGaussian(a);
			return value.GetAverageLog(value) - copy.GetAverageLog(copy);
		}

		/// <summary>
		/// VMP message to 'copy'
		/// </summary>
		/// <param name="value">Incoming message from 'value'. Must be a proper distribution.  If uniform, the result will be uniform.</param>
		/// <param name="to_value">Previous outgoing message to 'value'.</param>
		/// <returns>The outgoing VMP message to the 'copy' argument</returns>
		/// <remarks><para>
		/// The outgoing message is a distribution matching the moments of 'copy' as the random arguments are varied.
		/// The formula is <c>proj[sum_(value) p(value) factor(copy,value)]</c>.
		/// </para></remarks>
		/// <exception cref="ImproperMessageException"><paramref name="value"/> is not a proper distribution</exception>
		[Quality(QualityBand.Preview)]
		public static TruncatedGaussian CopyAverageLogarithm([SkipIfUniform, Stochastic] Gaussian value, Gaussian to_value)
		{
			return new TruncatedGaussian(value / to_value);
		}

		/*-------------- Nonconjugate Gaussian --------------------*/

		/// <summary>
		/// VMP message to 'value'
		/// </summary>
		/// <param name="copy">Incoming message from 'copy'.</param>
		/// <param name="value">Incoming message from 'value'.</param>
		/// <param name="to_value">Previous outgoing message to 'value'.</param>
		/// <returns>The outgoing VMP message to the 'value' argument</returns>
		/// <remarks><para>
		/// The outgoing message is the factor viewed as a function of 'value' with 'copy' integrated out.
		/// The formula is <c>sum_copy p(copy) factor(copy,value)</c>.
		/// </para><para>We reverse the direction of factor to get the behaviour we want here</para></remarks>
		public static Gaussian ValueAverageLogarithm(NonconjugateGaussian copy, Gaussian value, Gaussian to_value)
		{
			var a = value / to_value;
			copy *= new NonconjugateGaussian(a);
			var result = copy.GetGaussian(true) / a;
			return result;
		}

		/// <summary>
		/// VMP message to 'copy'
		/// </summary>
		/// <param name="value">Incoming message from 'value'.</param>
		/// <param name="copy">Incoming message from 'copy'.</param>
		/// <param name="to_value">Previous outgoing message to 'value'.</param>
		/// <returns>The outgoing VMP message to the 'copy' argument</returns>
		/// <remarks><para>
		/// The outgoing message is a distribution matching the moments of 'copy' as the random arguments are varied.
		/// The formula is <c>proj[sum_(value) p(value) factor(copy,value)]</c>.
		/// </para><para>This message should include the previous contribution.</para></remarks>
		public static NonconjugateGaussian CopyAverageLogarithm(Gaussian value, NonconjugateGaussian copy, Gaussian to_value)
		{
			return new NonconjugateGaussian(value / to_value);
		}
	}

	[FactorMethod(new string[] { "copy", "Value" }, typeof(Factor), "Copy<>")]
	[Quality(QualityBand.Experimental)]
	public class MaxProductCopyOp<T>
	{
		/// <summary>
		/// Max product message to 'value'
		/// </summary>
		/// <param name="copy">Constant value for 'copy'.</param>
		/// <param name="result">Modified to contain the outgoing message</param>
		/// <returns><paramref name="result"/></returns>
		/// <remarks><para>
		/// The outgoing message is the factor viewed as a function of 'value' conditioned on the given values.
		/// </para></remarks>
		public static TDist ValueMaxConditional<TDist>([IsReturned] TDist copy, TDist result)
			where TDist : IDistribution<T>
		{
			return copy;
		}

		/// <summary>
		/// Max product message to 'copy'
		/// </summary>
		/// <param name="Value">Incoming message from 'value'.</param>
		/// <returns>The outgoing VMP message to the 'copy' argument</returns>
		/// <remarks><para>
		/// The outgoing message is a distribution matching the moments of 'copy' as the random arguments are varied.
		/// The formula is <c>proj[sum_(value) p(value) factor(copy,value)]</c>.
		/// </para></remarks>
		public static TDist CopyMaxConditional<TDist>([IsReturned] TDist Value)
			where TDist : IDistribution<T>
		{
			return Value;
		}
	}
}
