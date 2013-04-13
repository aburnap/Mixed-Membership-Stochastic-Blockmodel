// (C) Copyright 2008 Microsoft Research Cambridge
using System;
using System.Collections.Generic;
using System.Text;
using MicrosoftResearch.Infer.Distributions;
using MicrosoftResearch.Infer.Maths;
using MicrosoftResearch.Infer.Utils;
using MicrosoftResearch.Infer.Collections;

namespace MicrosoftResearch.Infer.Factors
{
	/// <summary>
	/// Provides outgoing messages for <see cref="Factor.GetItem2D{T}"/>, given random arguments to the function.
	/// This factor gets an item from an array of items
	/// </summary>
	[FactorMethod(typeof(Factor), "GetItem2D<>")]
	[Quality(QualityBand.Experimental)]
	public static class GetItem2DOp<T>
	{
		/// <summary>
		/// Evidence message for EP
		/// </summary>
		/// <param name="item">Incoming message from 'item'.</param>
		/// <param name="to_item">Outgoing message to 'item'.</param>
		/// <returns>Logarithm of the factor's average value across the given argument distributions</returns>
		/// <remarks><para>
		/// The formula for the result is <c>log(sum_(item) p(item) factor(item,array,index1,index2))</c>.
		/// </para></remarks>
		public static double LogAverageFactor<Distribution>(Distribution item, [Fresh] Distribution to_item)
			where Distribution : CanGetLogAverageOf<Distribution>
		{
			return to_item.GetLogAverageOf(item);
		}
		/// <summary>
		/// Evidence message for EP
		/// </summary>
		/// <param name="item">Constant value for 'item'.</param>
		/// <param name="array">Incoming message from 'array'.</param>
		/// <param name="index1">Constant value for 'index1'.</param>
		/// <param name="index2">Constant value for 'index2'.</param>
		/// <returns>Logarithm of the factor's average value across the given argument distributions</returns>
		/// <remarks><para>
		/// The formula for the result is <c>log(sum_(array) p(array) factor(item,array,index1,index2))</c>.
		/// </para></remarks>
		public static double LogAverageFactor<Distribution>(T item, IArray2D<Distribution> array, int index1, int index2)
			where Distribution : CanGetLogProb<T>
		{
			return array[index1,index2].GetLogProb(item);
		}

		/// <summary>
		/// Evidence message for EP
		/// </summary>
		/// <param name="item">Incoming message from 'item'.</param>
		/// <returns>Logarithm of the factor's contribution the EP model evidence</returns>
		/// <remarks><para>
		/// The formula for the result is <c>log(sum_(item) p(item) factor(item,array,index1,index2) / sum_item p(item) messageTo(item))</c>.
		/// Adding up these values across all factors and variables gives the log-evidence estimate for EP.
		/// </para></remarks>
		[Skip]
		public static double LogEvidenceRatio<Distribution>(Distribution item) { return 0.0; }
		/// <summary>
		/// Evidence message for EP
		/// </summary>
		/// <param name="item">Constant value for 'item'.</param>
		/// <param name="array">Incoming message from 'array'.</param>
		/// <param name="index1">Constant value for 'index1'.</param>
		/// <param name="index2">Constant value for 'index2'.</param>
		/// <returns>Logarithm of the factor's contribution the EP model evidence</returns>
		/// <remarks><para>
		/// The formula for the result is <c>log(sum_(array) p(array) factor(item,array,index1,index2))</c>.
		/// Adding up these values across all factors and variables gives the log-evidence estimate for EP.
		/// </para></remarks>
		public static double LogEvidenceRatio<Distribution>(T item, IArray2D<Distribution> array, int index1, int index2)
			where Distribution : CanGetLogProb<T>
		{ return LogAverageFactor(item, array, index1, index2); }

		/// <summary>
		/// EP message to 'item'
		/// </summary>
		/// <param name="array">Incoming message from 'array'. Must be a proper distribution.  If all elements are uniform, the result will be uniform.</param>
		/// <param name="index1">Constant value for 'index1'.</param>
		/// <param name="index2">Constant value for 'index2'.</param>
		/// <param name="result">Modified to contain the outgoing message</param>
		/// <returns><paramref name="result"/></returns>
		/// <remarks><para>
		/// The outgoing message is a distribution matching the moments of 'item' as the random arguments are varied.
		/// The formula is <c>proj[p(item) sum_(array) p(array) factor(item,array,index1,index2)]/p(item)</c>.
		/// </para></remarks>
		/// <exception cref="ImproperMessageException"><paramref name="array"/> is not a proper distribution</exception>
		public static Distribution ItemAverageConditional<Distribution>([SkipIfAllUniform] IArray2D<Distribution> array, int index1, int index2, Distribution result)
			where Distribution : SettableTo<Distribution>
		{
			result.SetTo(array[index1,index2]);
			return result;
		}

		/// <summary>
		/// EP message to 'array'
		/// </summary>
		/// <param name="item">Incoming message from 'item'. Must be a proper distribution.  If uniform, the result will be uniform.</param>
		/// <param name="index1">Constant value for 'index1'.</param>
		/// <param name="index2">Constant value for 'index2'.</param>
		/// <param name="result">Modified to contain the outgoing message</param>
		/// <returns><paramref name="result"/></returns>
		/// <remarks><para>
		/// The outgoing message is a distribution matching the moments of 'array' as the random arguments are varied.
		/// The formula is <c>proj[p(array) sum_(item) p(item) factor(item,array,index1,index2)]/p(array)</c>.
		/// </para></remarks>
		/// <exception cref="ImproperMessageException"><paramref name="item"/> is not a proper distribution</exception>
		public static DistributionArray ArrayAverageConditional<Distribution, DistributionArray>([SkipIfUniform] Distribution item, int index1, int index2, DistributionArray result)
			where DistributionArray : IArray2D<Distribution>
			where Distribution : SettableTo<Distribution>
		{
			// assume result is initialized to uniform.
			Distribution value = result[index1,index2];
			value.SetTo(item);
			result[index1,index2] = value;
			return result;
		}

		/// <summary>
		/// EP message to 'array'
		/// </summary>
		/// <param name="item">Constant value for 'item'.</param>
		/// <param name="index1">Constant value for 'index1'.</param>
		/// <param name="index2">Constant value for 'index2'.</param>
		/// <param name="result">Modified to contain the outgoing message</param>
		/// <returns><paramref name="result"/></returns>
		/// <remarks><para>
		/// The outgoing message is the factor viewed as a function of 'array' conditioned on the given values.
		/// </para></remarks>
		public static DistributionArray ArrayAverageConditional<Distribution, DistributionArray>(T item, int index1, int index2, DistributionArray result)
			where DistributionArray : IArray2D<Distribution>
			where Distribution : HasPoint<T>
		{
			// assume result is initialized to uniform.
			Distribution value = result[index1,index2];
			value.Point = item;
			result[index1,index2] = value;
			return result;
		}

		//-- VMP -------------------------------------------------------------------------------------------------------------

		/// <summary>
		/// Evidence message for VMP
		/// </summary>
		/// <returns>Zero</returns>
		/// <remarks><para>
		/// The formula for the result is <c>log(factor(item,array,index1,index2))</c>.
		/// Adding up these values across all factors and variables gives the log-evidence estimate for VMP.
		/// </para></remarks>
		[Skip]
		public static double AverageLogFactor() { return 0.0; }

		/// <summary>
		/// VMP message to 'item'
		/// </summary>
		/// <param name="array">Incoming message from 'array'. Must be a proper distribution.  If all elements are uniform, the result will be uniform.</param>
		/// <param name="index1">Constant value for 'index1'.</param>
		/// <param name="index2">Constant value for 'index2'.</param>
		/// <param name="result">Modified to contain the outgoing message</param>
		/// <returns><paramref name="result"/></returns>
		/// <remarks><para>
		/// The outgoing message is a distribution matching the moments of 'item' as the random arguments are varied.
		/// The formula is <c>proj[sum_(array) p(array) factor(item,array,index1,index2)]</c>.
		/// </para></remarks>
		/// <exception cref="ImproperMessageException"><paramref name="array"/> is not a proper distribution</exception>
		public static Distribution ItemAverageLogarithm<Distribution>([SkipIfAllUniform] IArray2D<Distribution> array, int index1, int index2, Distribution result)
			where Distribution : SettableTo<Distribution>
		{
			result.SetTo(array[index1,index2]);
			return result;
		}

		/// <summary>
		/// VMP message to 'array'
		/// </summary>
		/// <param name="item">Incoming message from 'item'. Must be a proper distribution.  If uniform, the result will be uniform.</param>
		/// <param name="index1">Constant value for 'index1'.</param>
		/// <param name="index2">Constant value for 'index2'.</param>
		/// <param name="result">Modified to contain the outgoing message</param>
		/// <returns><paramref name="result"/></returns>
		/// <remarks><para>
		/// The outgoing message is the factor viewed as a function of 'array' with 'item' integrated out.
		/// The formula is <c>sum_item p(item) factor(item,array,index1,index2)</c>.
		/// </para></remarks>
		/// <exception cref="ImproperMessageException"><paramref name="item"/> is not a proper distribution</exception>
		public static DistributionArray ArrayAverageLogarithm<Distribution, DistributionArray>([SkipIfUniform] Distribution item, int index1, int index2, DistributionArray result)
			where DistributionArray : IArray2D<Distribution>
			where Distribution : SettableTo<Distribution>
		{
			Distribution value = result[index1,index2];
			value.SetTo(item);
			result[index1,index2] = value;
			return result;
		}

		/// <summary>
		/// VMP message to 'array'
		/// </summary>
		/// <param name="item">Constant value for 'item'.</param>
		/// <param name="index1">Constant value for 'index1'.</param>
		/// <param name="index2">Constant value for 'index2'.</param>
		/// <param name="result">Modified to contain the outgoing message</param>
		/// <returns><paramref name="result"/></returns>
		/// <remarks><para>
		/// The outgoing message is the factor viewed as a function of 'array' conditioned on the given values.
		/// </para></remarks>
		public static DistributionArray ArrayAverageLogarithm<Distribution, DistributionArray>(T item, int index1, int index2, DistributionArray result)
			where DistributionArray : IArray2D<Distribution>
			where Distribution : HasPoint<T>
		{
			// assume result is initialized to uniform.
			Distribution value = result[index1,index2];
			value.Point = item;
			result[index1,index2] = value;
			return result;
		}
	}

}
