// (C) Copyright 2008 Microsoft Research Cambridge
using System;
using System.Collections.Generic;
using System.Text;
using MicrosoftResearch.Infer.Distributions;
using MicrosoftResearch.Infer.Maths;
using MicrosoftResearch.Infer.Utils;

namespace MicrosoftResearch.Infer.Factors
{
	/// <summary>
	/// Provides outgoing messages for <see cref="Factor.GetItem{T}"/>, given random arguments to the function.
	/// This factor gets an item from an array of items
	/// </summary>
	[FactorMethod(typeof(Factor), "GetItem<>")]
	[Quality(QualityBand.Preview)]
	public static class GetItemOp<T>
	{
		public static double LogAverageFactor(T item, IList<T> array, int index)
		{
			IEqualityComparer<T> equalityComparer = Utils.Util.GetEqualityComparer<T>();
			return equalityComparer.Equals(item,array[index]) ? 0.0 : Double.NegativeInfinity;
		}
		public static double LogEvidenceRatio(T item, IList<T> array, int index) { return LogAverageFactor(item, array, index); }
		public static double AverageLogFactor(T item, IList<T> array, int index) { return LogAverageFactor(item, array, index); }

		/// <summary>
		/// Evidence message for EP
		/// </summary>
		/// <param name="item">Incoming message from 'item'.</param>
		/// <param name="to_item">Outgoing message to 'item'.</param>
		/// <returns>Logarithm of the factor's average value across the given argument distributions</returns>
		/// <remarks><para>
		/// The formula for the result is <c>log(sum_(item) p(item) factor(item,array,index))</c>.
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
		/// <param name="index">Constant value for 'index'.</param>
		/// <returns>Logarithm of the factor's average value across the given argument distributions</returns>
		/// <remarks><para>
		/// The formula for the result is <c>log(sum_(array) p(array) factor(item,array,index))</c>.
		/// </para></remarks>
		public static double LogAverageFactor<Distribution>(T item, IList<Distribution> array, int index)
			where Distribution : CanGetLogProb<T>
		{
			return array[index].GetLogProb(item);
		}

		/// <summary>
		/// Evidence message for EP
		/// </summary>
		/// <param name="item">Incoming message from 'item'.</param>
		/// <returns>Logarithm of the factor's contribution the EP model evidence</returns>
		/// <remarks><para>
		/// The formula for the result is <c>log(sum_(item) p(item) factor(item,array,index) / sum_item p(item) messageTo(item))</c>.
		/// Adding up these values across all factors and variables gives the log-evidence estimate for EP.
		/// </para></remarks>
		[Skip]
		public static double LogEvidenceRatio<Distribution>(Distribution item) where Distribution : IDistribution<T>
		{ return 0.0; }
		/// <summary>
		/// Evidence message for EP
		/// </summary>
		/// <param name="item">Constant value for 'item'.</param>
		/// <param name="array">Incoming message from 'array'.</param>
		/// <param name="index">Constant value for 'index'.</param>
		/// <returns>Logarithm of the factor's contribution the EP model evidence</returns>
		/// <remarks><para>
		/// The formula for the result is <c>log(sum_(array) p(array) factor(item,array,index))</c>.
		/// Adding up these values across all factors and variables gives the log-evidence estimate for EP.
		/// </para></remarks>
		public static double LogEvidenceRatio<Distribution>(T item, IList<Distribution> array, int index)
			where Distribution : CanGetLogProb<T>
		{ return LogAverageFactor(item, array, index); }

		/// <summary>
		/// EP message to 'item'
		/// </summary>
		/// <param name="array">Incoming message from 'array'. Must be a proper distribution.  If all elements are uniform, the result will be uniform.</param>
		/// <param name="index">Constant value for 'index'.</param>
		/// <param name="result">Modified to contain the outgoing message</param>
		/// <returns><paramref name="result"/></returns>
		/// <remarks><para>
		/// The outgoing message is a distribution matching the moments of 'item' as the random arguments are varied.
		/// The formula is <c>proj[p(item) sum_(array) p(array) factor(item,array,index)]/p(item)</c>.
		/// </para></remarks>
		/// <exception cref="ImproperMessageException"><paramref name="array"/> is not a proper distribution</exception>
		public static Distribution ItemAverageConditional<Distribution>([SkipIfAllUniform] IList<Distribution> array, int index, Distribution result)
			where Distribution : SettableTo<Distribution>
		{
			result.SetTo(array[index]);
			return result;
		}
		[Skip]
		public static Distribution ItemAverageConditionalInit<Distribution>([IgnoreDependency] IList<Distribution> array)
			where Distribution : ICloneable
		{
			return (Distribution)array[0].Clone();
		}

		/// <summary>
		/// EP message to 'array'
		/// </summary>
		/// <param name="item">Incoming message from 'item'. Must be a proper distribution.  If uniform, the result will be uniform.</param>
		/// <param name="index">Constant value for 'index'.</param>
		/// <param name="result">Modified to contain the outgoing message</param>
		/// <returns><paramref name="result"/></returns>
		/// <remarks><para>
		/// The outgoing message is a distribution matching the moments of 'array' as the random arguments are varied.
		/// The formula is <c>proj[p(array) sum_(item) p(item) factor(item,array,index)]/p(array)</c>.
		/// </para></remarks>
		/// <exception cref="ImproperMessageException"><paramref name="item"/> is not a proper distribution</exception>
		public static DistributionArray ArrayAverageConditional<Distribution, DistributionArray>([SkipIfUniform] Distribution item, int index, DistributionArray result)
			where DistributionArray : IList<Distribution>
			where Distribution : SettableTo<Distribution>
		{
			// assume result is initialized to uniform.
			Distribution value = result[index];
			value.SetTo(item);
			result[index] = value;
			return result;
		}

		/// <summary>
		/// EP message to 'array'
		/// </summary>
		/// <param name="item">Constant value for 'item'.</param>
		/// <param name="index">Constant value for 'index'.</param>
		/// <param name="result">Modified to contain the outgoing message</param>
		/// <returns><paramref name="result"/></returns>
		/// <remarks><para>
		/// The outgoing message is the factor viewed as a function of 'array' conditioned on the given values.
		/// </para></remarks>
		public static DistributionArray ArrayAverageConditional<Distribution, DistributionArray>(T item, int index, DistributionArray result)
			where DistributionArray : IList<Distribution>
			where Distribution : HasPoint<T>
		{
			// assume result is initialized to uniform.
			Distribution value = result[index];
			value.Point = item;
			result[index] = value;
			return result;
		}

		//-- VMP -------------------------------------------------------------------------------------------------------------

		/// <summary>
		/// Evidence message for VMP
		/// </summary>
		/// <returns>Zero</returns>
		/// <remarks><para>
		/// The formula for the result is <c>log(factor(item,array,index))</c>.
		/// Adding up these values across all factors and variables gives the log-evidence estimate for VMP.
		/// </para></remarks>
		[Skip]
		public static double AverageLogFactor() { return 0.0; }

		/// <summary>
		/// VMP message to 'item'
		/// </summary>
		/// <param name="array">Incoming message from 'array'. Must be a proper distribution.  If all elements are uniform, the result will be uniform.</param>
		/// <param name="index">Constant value for 'index'.</param>
		/// <param name="result">Modified to contain the outgoing message</param>
		/// <returns><paramref name="result"/></returns>
		/// <remarks><para>
		/// The outgoing message is a distribution matching the moments of 'item' as the random arguments are varied.
		/// The formula is <c>proj[sum_(array) p(array) factor(item,array,index)]</c>.
		/// </para></remarks>
		/// <exception cref="ImproperMessageException"><paramref name="array"/> is not a proper distribution</exception>
		public static Distribution ItemAverageLogarithm<Distribution>([SkipIfAllUniform] IList<Distribution> array, int index, Distribution result)
			where Distribution : SettableTo<Distribution>
		{
			result.SetTo(array[index]);
			return result;
		}
		[Skip]
		public static Distribution ItemAverageLogarithmInit<Distribution>([IgnoreDependency] IList<Distribution> array)
			where Distribution : ICloneable
		{
			return (Distribution)array[0].Clone();
		}

		/// <summary>
		/// VMP message to 'array'
		/// </summary>
		/// <param name="item">Incoming message from 'item'. Must be a proper distribution.  If uniform, the result will be uniform.</param>
		/// <param name="index">Constant value for 'index'.</param>
		/// <param name="result">Modified to contain the outgoing message</param>
		/// <returns><paramref name="result"/></returns>
		/// <remarks><para>
		/// The outgoing message is the factor viewed as a function of 'array' with 'item' integrated out.
		/// The formula is <c>sum_item p(item) factor(item,array,index)</c>.
		/// </para></remarks>
		/// <exception cref="ImproperMessageException"><paramref name="item"/> is not a proper distribution</exception>
		public static DistributionArray ArrayAverageLogarithm<Distribution, DistributionArray>([SkipIfUniform] Distribution item, int index, DistributionArray result)
			where DistributionArray : IList<Distribution>
			where Distribution : SettableTo<Distribution>
		{
			Distribution value = result[index];
			value.SetTo(item);
			result[index] = value;
			return result;
		}

		/// <summary>
		/// VMP message to 'array'
		/// </summary>
		/// <param name="item">Constant value for 'item'.</param>
		/// <param name="index">Constant value for 'index'.</param>
		/// <param name="result">Modified to contain the outgoing message</param>
		/// <returns><paramref name="result"/></returns>
		/// <remarks><para>
		/// The outgoing message is the factor viewed as a function of 'array' conditioned on the given values.
		/// </para></remarks>
		public static DistributionArray ArrayAverageLogarithm<Distribution, DistributionArray>(T item, int index, DistributionArray result)
			where DistributionArray : IList<Distribution>
			where Distribution : HasPoint<T>
		{
			// assume result is initialized to uniform.
			Distribution value = result[index];
			value.Point = item;
			result[index] = value;
			return result;
		}
	}

}
