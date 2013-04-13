using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using MicrosoftResearch.Infer.Distributions;
using MicrosoftResearch.Infer.Maths;
using MicrosoftResearch.Infer.Utils;

namespace MicrosoftResearch.Infer.Factors
{
	[FactorMethod(typeof(Factor), "JaggedSubarray<>")]
	[Buffers("marginal")]
	[Quality(QualityBand.Preview)]
	public static class JaggedSubarrayOp<T>
	{
		/// <summary>
		/// Evidence message for EP
		/// </summary>
		/// <param name="items">Incoming message from 'items'.</param>
		/// <param name="array">Constant value for 'array'.</param>
		/// <param name="indices">Constant value for 'indices'.</param>
		/// <returns>Logarithm of the factor's average value across the given argument distributions</returns>
		/// <remarks><para>
		/// The formula for the result is <c>log(sum_(items) p(items) factor(items,array,indices))</c>.
		/// </para></remarks>
		public static double LogAverageFactor<ItemType>(IList<ItemType> items, IList<T> array, IList<int[]> indices)
			where ItemType : IList<T>
		{
			IEqualityComparer<T> equalityComparer = Utils.Util.GetEqualityComparer<T>();
			for (int i = 0; i < indices.Count; i++) {
				for (int j = 0; j < indices[i].Length; j++) {
					if (!equalityComparer.Equals(items[i][j],array[indices[i][j]])) return Double.NegativeInfinity;
				}
			}
			return 0.0;
		}
		public static double LogEvidenceRatio<ItemType>(IList<ItemType> items, IList<T> array, IList<int[]> indices)
			where ItemType : IList<T>
		{ return LogAverageFactor<ItemType>(items, array, indices); }
		public static double AverageLogFactor<ItemType>(IList<ItemType> items, IList<T> array, IList<int[]> indices)
			where ItemType : IList<T>
		{ return LogAverageFactor<ItemType>(items, array, indices); }

		/// <summary>
		/// Evidence message for EP
		/// </summary>
		/// <param name="items">Incoming message from 'items'.</param>
		/// <param name="array">Incoming message from 'array'.</param>
		/// <param name="indices">Constant value for 'indices'.</param>
		/// <returns>Logarithm of the factor's average value across the given argument distributions</returns>
		/// <remarks><para>
		/// The formula for the result is <c>log(sum_(items,array) p(items,array) factor(items,array,indices))</c>.
		/// </para></remarks>
		public static double LogAverageFactor<DistributionType, ItemType>(IList<ItemType> items, IList<DistributionType> array, IList<int[]> indices)
			where DistributionType : IDistribution<T>, SettableToProduct<DistributionType>, CanGetLogAverageOf<DistributionType>
			where ItemType : IList<DistributionType>
		{
			double z = 0.0;
			Dictionary<int,DistributionType> productBefore = new Dictionary<int, DistributionType>();
			for (int i = 0; i < indices.Count; i++) {
				for (int j = 0; j < indices[i].Length; j++) {
					DistributionType value;
					if (!productBefore.TryGetValue(indices[i][j], out value)) {
						value = (DistributionType)array[indices[i][j]].Clone();
					}
					z += value.GetLogAverageOf(items[i][j]);
					value.SetToProduct(value, items[i][j]);
					productBefore[indices[i][j]] = value;
				}
			}
			return z;
		}
		[Skip]
		public static double AverageLogFactor<DistributionType, ItemType>(IList<ItemType> items, IList<DistributionType> array, IList<int[]> indices)
			where DistributionType : IDistribution<T>, SettableToProduct<DistributionType>, CanGetLogAverageOf<DistributionType>
			where ItemType : IList<DistributionType>
		{
			return 0.0;
		}

		/// <summary>
		/// Evidence message for EP
		/// </summary>
		/// <param name="items">Incoming message from 'items'.</param>
		/// <param name="array">Incoming message from 'array'.</param>
		/// <param name="indices">Constant value for 'indices'.</param>
		/// <returns>Logarithm of the factor's average value across the given argument distributions</returns>
		/// <remarks><para>
		/// The formula for the result is <c>log(sum_(items,array) p(items,array) factor(items,array,indices))</c>.
		/// </para></remarks>
		public static double LogAverageFactor<DistributionType, ItemType>(IList<DistributionType> array, IList<ItemType> items, IList<int[]> indices)
			where DistributionType : IDistribution<T>, CanGetLogProb<T>
			where ItemType : IList<T>
		{
			double z = 0.0;
			Dictionary<int,DistributionType> productBefore = new Dictionary<int, DistributionType>();
			for (int i = 0; i < indices.Count; i++) {
				for (int j = 0; j < indices[i].Length; j++) {
					DistributionType value;
					if(!productBefore.TryGetValue(indices[i][j], out value)) {
						value = array[indices[i][j]];
					}
					z += value.GetLogProb(items[i][j]);
					value.Point = items[i][j];
					productBefore[indices[i][j]] = value;
				}
			}
			return z;
		}
		/// <summary>
		/// Evidence message for EP
		/// </summary>
		/// <param name="items">Incoming message from 'items'.</param>
		/// <param name="array">Incoming message from 'array'.</param>
		/// <param name="indices">Constant value for 'indices'.</param>
		/// <returns>Logarithm of the factor's contribution the EP model evidence</returns>
		/// <remarks><para>
		/// The formula for the result is <c>log(sum_(items,array) p(items,array) factor(items,array,indices) / sum_items p(items) messageTo(items))</c>.
		/// Adding up these values across all factors and variables gives the log-evidence estimate for EP.
		/// </para></remarks>
		public static double LogEvidenceRatio<DistributionType, ItemType>(IList<DistributionType> array, IList<ItemType> items, IList<int[]> indices)
			where DistributionType : IDistribution<T>, CanGetLogProb<T>
			where ItemType : IList<T>
		{
			return LogAverageFactor<DistributionType, ItemType>(array, items, indices);
		}
		[Skip]
		public static double AverageLogFactor<DistributionType, ItemType>(IList<DistributionType> array, IList<ItemType> items, IList<int[]> indices)
			where DistributionType : IDistribution<T>, CanGetLogProb<T>
			where ItemType : IList<T>
		{
			return 0.0;
		}

		/// <summary>
		/// Evidence message for EP
		/// </summary>
		/// <param name="items">Incoming message from 'items'.</param>
		/// <param name="array">Constant value for 'array'.</param>
		/// <param name="indices">Constant value for 'indices'.</param>
		/// <returns>Logarithm of the factor's average value across the given argument distributions</returns>
		/// <remarks><para>
		/// The formula for the result is <c>log(sum_(items) p(items) factor(items,array,indices))</c>.
		/// </para></remarks>
		public static double LogAverageFactor<DistributionType, ItemType>(IList<ItemType> items, IList<T> array, IList<int[]> indices)
			where DistributionType : CanGetLogProb<T>
			where ItemType : IList<DistributionType>
		{
			double z = 0.0;
			for (int i = 0; i < indices.Count; i++)
				for (int j = 0; j < indices[i].Length; j++)
					z += items[i][j].GetLogProb(array[indices[i][j]]);
			return z;
		}
		[Skip]
		public static double LogEvidenceRatio<DistributionType, ItemType>(IList<ItemType> items, IList<T> array, IList<int[]> indices)
			where DistributionType : CanGetLogProb<T>
			where ItemType : IList<DistributionType>
		{
			return 0.0;
		}
		[Skip]
		public static double AverageLogFactor<DistributionType, ItemType>(IList<ItemType> items, IList<T> array, IList<int[]> indices)
			where DistributionType : CanGetLogProb<T>
			where ItemType : IList<DistributionType>
		{
			return 0.0;
		}

#if true
		/// <summary>
		/// Evidence message for EP
		/// </summary>
		/// <param name="items">Incoming message from 'items'.</param>
		/// <param name="array">Incoming message from 'array'.</param>
		/// <param name="indices">Constant value for 'indices'.</param>
		/// <returns>Logarithm of the factor's contribution the EP model evidence</returns>
		/// <remarks><para>
		/// The formula for the result is <c>log(sum_(items,array) p(items,array) factor(items,array,indices) / sum_items p(items) messageTo(items))</c>.
		/// Adding up these values across all factors and variables gives the log-evidence estimate for EP.
		/// </para></remarks>
		public static double LogEvidenceRatio<DistributionType, ItemArrayType, ItemType>(ItemArrayType items, IList<DistributionType> array, IList<int[]> indices, ItemArrayType to_items)
			where ItemArrayType : IList<ItemType>, ICloneable
			where ItemType : IList<DistributionType>
			where DistributionType : SettableToUniform, SettableToProduct<DistributionType>, CanGetLogAverageOf<DistributionType>, ICloneable
		{
			// this code is adapted from GetItemsOp
			double z = 0.0;
			if (items.Count == 0) return 0.0;
			if (items.Count == 1 && items[0].Count <= 1) return 0.0;
			Dictionary<int,DistributionType> productBefore = new Dictionary<int, DistributionType>();
			for (int i = 0; i < indices.Count; i++) {
				for (int j = 0; j < indices[i].Length; j++) {
					DistributionType value;
					if (!productBefore.TryGetValue(indices[i][j], out value)) {
						value = (DistributionType)array[indices[i][j]].Clone();
					}
					z += value.GetLogAverageOf(items[i][j]);
					value.SetToProduct(value, items[i][j]);
					productBefore[indices[i][j]] = value;
					z -= to_items[i][j].GetLogAverageOf(items[i][j]);
				}
			}
			return z;
		}
#else
		/// <summary>
		/// Evidence message for EP
		/// </summary>
		/// <param name="items">Incoming message from 'items'.</param>
		/// <param name="array">Incoming message from 'array'.</param>
		/// <param name="indices">Constant value for 'indices'.</param>
		/// <returns>Logarithm of the factor's contribution the EP model evidence</returns>
		/// <remarks><para>
		/// The formula for the result is <c>log(sum_(items,array) p(items,array) factor(items,array,indices) / sum_items p(items) messageTo(items))</c>.
		/// Adding up these values across all factors and variables gives the log-evidence estimate for EP.
		/// </para></remarks>
		public static double LogEvidenceRatio<DistributionType, ItemArrayType, ItemType>(ItemArrayType items, IList<DistributionType> array, IList<int[]> indices)
			where ItemArrayType : IList<ItemType>, ICloneable
			where ItemType : IList<DistributionType>
			where DistributionType : SettableToUniform, SettableToProduct<DistributionType>, CanGetLogAverageOf<DistributionType>, ICloneable
		{
			// this code is adapted from GetItemsOp
			double z = 0.0;
			if (items.Count == 0) return 0.0;
			if (items.Count == 1 && items[0].Count <= 1) return 0.0;
			bool firstTime = true;
			DistributionType uniform = default(DistributionType);
			ItemArrayType productBefore = (ItemArrayType)items.Clone();
			ItemArrayType productAfter = (ItemArrayType)items.Clone();
			Dictionary<int,KeyValuePair<int,int>> indexToItem = new Dictionary<int, KeyValuePair<int, int>>();
			for (int i = 0; i < indices.Count; i++) {
				for (int j = 0; j < indices[i].Length; j++) {
					KeyValuePair<int,int> previousItem;
					if (!indexToItem.TryGetValue(indices[i][j], out previousItem)) {
						// no previous item with this index
						productBefore[i][j] = array[indices[i][j]];
					} else {
						DistributionType temp = productBefore[i][j];
						temp.SetToProduct(productBefore[previousItem.Key][previousItem.Value], items[previousItem.Key][previousItem.Value]);
						productBefore[i][j] = temp;
					}
					z += productBefore[i][j].GetLogAverageOf(items[i][j]);
					indexToItem[indices[i][j]] = new KeyValuePair<int, int>(i, j);
					if (firstTime) {
						uniform = (DistributionType)items[i][j].Clone();
						uniform.SetToUniform();
					}
				}
			}
			indexToItem.Clear();
			for (int i = indices.Count - 1; i >= 0; i--) {
				for (int j = indices[i].Length - 1; j >= 0; j--) {
					KeyValuePair<int,int> itemAfter;
					if (!indexToItem.TryGetValue(indices[i][j], out itemAfter)) {
						// no item after with this index
						productAfter[i][j] = uniform;
					} else {
						DistributionType temp = productAfter[i][j];
						temp.SetToProduct(productAfter[itemAfter.Key][itemAfter.Value], items[itemAfter.Key][itemAfter.Value]);
						productAfter[i][j] = temp;
					}
					DistributionType toItem = (DistributionType)items[i][j].Clone();
					toItem.SetToProduct(productBefore[i][j], productAfter[i][j]);
					z -= toItem.GetLogAverageOf(items[i][j]);
					indexToItem[indices[i][j]] = new KeyValuePair<int,int>(i,j);
				}
			}
			return z;
		}
#endif

		[Skip]
		public static ArrayType MarginalInit<ArrayType>([SkipIfUniform] ArrayType array)
			where ArrayType : ICloneable
		{
			return (ArrayType)array.Clone();
		}
		public static ArrayType Marginal<ArrayType, DistributionType, ItemArrayType, ItemType>(ArrayType array, IList<ItemType> items, IList<int[]> indices, ArrayType result)
			where ArrayType : IList<DistributionType>, SettableTo<ArrayType>
			where DistributionType : SettableToProduct<DistributionType>
			where ItemType : IList<DistributionType>
		{
			Assert.IsTrue(items.Count == indices.Count, "items.Count != indices.Count");
			result.SetTo(array);
			for (int i = 0; i < indices.Count; i++) {
				for (int j = 0; j < indices[i].Length; j++) {
					DistributionType value = result[indices[i][j]];
					value.SetToProduct(value, items[i][j]);
					result[indices[i][j]] = value;
				}
			}
			return result;
		}
		public static ArrayType MarginalIncrement<ArrayType, DistributionType, ItemType>(ArrayType result, ItemType to_item, ItemType item, IList<int[]> indices, int resultIndex)
			where ArrayType : IList<DistributionType>, SettableTo<ArrayType>
			where DistributionType : SettableToProduct<DistributionType>
			where ItemType : IList<DistributionType>
		{
			int i = resultIndex;
			for (int j = 0; j < indices[i].Length; j++)	{
				DistributionType value = result[indices[i][j]];
				value.SetToProduct(to_item[j], item[j]);
				result[indices[i][j]] = value;
			}
			return result;
		}

		// items dependency must be ignored for Sequential schedule
		public static ItemType ItemsAverageConditional<ArrayType, DistributionType, ItemType>([Indexed, Cancels] ItemType items, [IgnoreDependency] ArrayType array, [SkipIfAllUniform] ArrayType marginal, IList<int[]> indices, int resultIndex, ItemType result)
			where ArrayType : IList<DistributionType>
			where ItemType : IList<DistributionType>
			where DistributionType : SettableToProduct<DistributionType>, SettableToRatio<DistributionType>
		{
			int i = resultIndex;
			Assert.IsTrue(result.Count == indices[i].Length, "result.Count != indices[i].Count");
			for (int j = 0; j < indices[i].Length; j++) {
				DistributionType value = result[j];
				value.SetToRatio(marginal[indices[i][j]], items[j]);
				result[j] = value;
			}
			return result;
		}

		public static ArrayType ArrayAverageConditional<DistributionType, ArrayType, ItemType>([SkipIfAllUniform] IList<ItemType> items, IList<int[]> indices, ArrayType result)
			where ArrayType : IList<DistributionType>, SettableToUniform
			where ItemType : IList<DistributionType>
			where DistributionType : SettableToProduct<DistributionType>
		{
			Assert.IsTrue(items.Count == indices.Count, "items.Count != indices.Count");
			result.SetToUniform();
			for (int i = 0; i < indices.Count; i++) {
				for (int j = 0; j < indices[i].Length; j++) {
					DistributionType value = result[indices[i][j]];
					value.SetToProduct(value, items[i][j]);
					result[indices[i][j]] = value;
				}
			}
			return result;
		}

		public static ArrayType ArrayAverageConditional<DistributionType, ArrayType, ItemType>(IList<int[]> indices, IList<ItemType> items, ArrayType result)
			where ArrayType : IList<DistributionType>, SettableToUniform
			where ItemType : IList<T>
			where DistributionType : HasPoint<T>
		{
			Assert.IsTrue(items.Count == indices.Count, "items.Count != indices.Count");
			result.SetToUniform();
			for (int i = 0; i < indices.Count; i++) {
				for (int j = 0; j < indices[i].Length; j++) {
					DistributionType value = result[indices[i][j]];
					value.Point = items[i][j];
					result[indices[i][j]] = value;
				}
			}
			return result;
		}

		//-- VMP -------------------------------------------------------------------------------------------------------------

		/// <summary>
		/// VMP message to 'items'
		/// </summary>
		/// <param name="array">Incoming message from 'array'. Must be a proper distribution.  If all elements are uniform, the result will be uniform.</param>
		/// <param name="indices">Constant value for 'indices'.</param>
		/// <param name="result">Modified to contain the outgoing message</param>
		/// <returns><paramref name="result"/></returns>
		/// <remarks><para>
		/// The outgoing message is a distribution matching the moments of 'items' as the random arguments are varied.
		/// The formula is <c>proj[sum_(array) p(array) factor(items,array,indices)]</c>.
		/// </para></remarks>
		/// <exception cref="ImproperMessageException"><paramref name="array"/> is not a proper distribution</exception>
		public static ResultType ItemsAverageLogarithm<DistributionType, ItemType, ResultType>([SkipIfAllUniform] IList<DistributionType> array, IList<int[]> indices, ResultType result)
			where ResultType : IList<ItemType>
			where ItemType : IList<DistributionType>
			where DistributionType : SettableTo<DistributionType>
		{
			for (int i = 0; i < indices.Count; i++) {
				for (int j = 0; j < indices[i].Length; j++) {
					DistributionType value = result[i][j];
					value.SetTo(array[indices[i][j]]);
					result[i][j] = value;
				}
			}
			return result;
		}

		/// <summary>
		/// VMP message to 'array'
		/// </summary>
		/// <param name="items">Incoming message from 'items'. Must be a proper distribution.  If all elements are uniform, the result will be uniform.</param>
		/// <param name="indices">Constant value for 'indices'.</param>
		/// <param name="result">Modified to contain the outgoing message</param>
		/// <returns><paramref name="result"/></returns>
		/// <remarks><para>
		/// The outgoing message is the factor viewed as a function of 'array' with 'items' integrated out.
		/// The formula is <c>sum_items p(items) factor(items,array,indices)</c>.
		/// </para></remarks>
		/// <exception cref="ImproperMessageException"><paramref name="items"/> is not a proper distribution</exception>
		public static ArrayType ArrayAverageLogarithm<DistributionType, ArrayType, ItemType>([SkipIfAllUniform] IList<ItemType> items, IList<int[]> indices, ArrayType result)
			where ArrayType : IList<DistributionType>, SettableToUniform
			where ItemType : IList<DistributionType>
			where DistributionType : SettableToUniform, SettableToProduct<DistributionType>
		{
			Assert.IsTrue(items.Count == indices.Count, "items.Count != indices.Count");
			result.SetToUniform();
			for (int i = 0; i < indices.Count; i++) {
				for (int j = 0; j < indices[i].Length; j++) {
					DistributionType value = result[indices[i][j]];
					value.SetToProduct(value, items[i][j]);
					result[indices[i][j]] = value;
				}
			}
			return result;
		}

		/// <summary>
		/// VMP message to 'array'
		/// </summary>
		/// <param name="items">Incoming message from 'items'.</param>
		/// <param name="indices">Constant value for 'indices'.</param>
		/// <param name="result">Modified to contain the outgoing message</param>
		/// <returns><paramref name="result"/></returns>
		/// <remarks><para>
		/// The outgoing message is the factor viewed as a function of 'array' with 'items' integrated out.
		/// The formula is <c>sum_items p(items) factor(items,array,indices)</c>.
		/// </para></remarks>
		public static ArrayType ArrayAverageLogarithm<DistributionType, ArrayType, ItemType>(IList<int[]> indices, IList<ItemType> items, ArrayType result)
			where ArrayType : IList<DistributionType>, SettableToUniform
			where ItemType : IList<T>
			where DistributionType : HasPoint<T>
		{
			return ArrayAverageConditional<DistributionType, ArrayType, ItemType>(indices, items, result);
		}
	}
}
