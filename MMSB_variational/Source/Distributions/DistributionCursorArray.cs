// (C) Copyright 2008 Microsoft Research Cambridge
using System;
using System.Collections.Generic;
using System.Text;
using MicrosoftResearch.Infer.Collections;
using MicrosoftResearch.Infer.Maths;
using MicrosoftResearch.Infer.Utils;

namespace MicrosoftResearch.Infer.Distributions
{
#if false
	/// <summary>
	/// An array of distribution objects, treated as values.
	/// </summary>
	/// <typeparam name="T">A distribution type.</typeparam>
	/// <remarks>
	/// Even if the distribution objects are reference types, they are treated as values during cloning
	/// and mutation.
	/// </remarks>
	public class DistributionArray<T> : ArrayAsList<T>, ICloneable,
		SettableTo<T>, SettableTo<DistributionArray<T>>, SettableToUniform,
		Diffable,
		SettableToProduct<DistributionArray<T>>,
		SettableToRatio<DistributionArray<T>>,
		SettableToPower<DistributionArray<T>>,
		SettableToWeightedSum<DistributionArray<T>>,
		CanGetLogAverageOf<DistributionArray<T>>,
		CanGetAverageLog<DistributionArray<T>>
		where T : ICloneable, SettableTo<T>, SettableToUniform,
		Diffable,
		SettableToProduct<T>,
		SettableToRatio<T>,
		SettableToPower<T>,
		SettableToWeightedSum<T>,
		CanGetLogAverageOf<T>,
		CanGetAverageLog<T>
	{
		public bool IsValueType;

#if false
		public static DistributionArray<T> CloneElementsOf(Array array)
		{
			DistributionArray<T> result = new DistributionArray<T>(StringUtil.ArrayDimensions(array));
			result.ModifyAll(array, delegate(T item, T aItem)
			{
				return (T)aItem.Clone();
			});
			return result;
		}
		public static DistributionArray<T> CloneElementsOfCollection(ICollection<T> array)
		{
			DistributionArray<T> result = new DistributionArray<T>(array.Count);
			int i = 0;
			foreach (T item in array) {
				result.array1D[i++] = (T)item.Clone();
			}
			return result;
		}
		public static DistributionArray<T> FromLength(params int[] lengths)
		{
			Type domainType = Distribution.GetDomainType(typeof(T));
			Type arrayType = Util.MakeArrayType(typeof(T), lengths.Length);
			Type type = typeof(DistributionArray<,,>).MakeGenericType(arrayType, typeof(T), domainType);
			return (DistributionArray<T>)Activator.CreateInstance(type, lengths);
		}
		public static DistributionArray<T> FromCloning(T value, params int[] lengths)
		{
			DistributionArray<T> result = FromLength(lengths);
			result.InitializeTo(value);
			return result;
		}
#endif

		public DistributionArray(Array array)
			: base(array)
		{
			IsValueType = typeof(T).IsValueType;
		}
		public DistributionArray(int length)
			: base(length)
		{
			IsValueType = typeof(T).IsValueType;
		}
		public DistributionArray(int length0, int length1)
			: base(length0, length1)
		{
			IsValueType = typeof(T).IsValueType;
		}
		public DistributionArray(params int[] lengths)
			: base(lengths)
		{
			IsValueType = typeof(T).IsValueType;
		}

		public DistributionArray(T value, int length)
			: this(length)
		{
			InitializeTo(value);
		}
		public DistributionArray(T value, params int[] lengths)
			: this(lengths)
		{
			InitializeTo(value);
		}
		/// <summary>
		/// Copy constructor.
		/// </summary>
		/// <param name="that"></param>
		public DistributionArray(DistributionArray<T> that)
			: this(that.GetLengths())
		{
			InitializeTo(that);
		}

		public void InitializeTo(T value)
		{
			if (IsValueType) {
				//ModifyAll(delegate(T item) { return value; });
				if (array1D != null) {
					for (int i = 0; i < array1D.Length; i++) {
						array1D[i] = value;
					}
				} else if (array2D != null) {
					for (int i = 0; i < array2D.GetLength(0); i++) {
						for (int j = 0; j < array2D.GetLength(1); j++) {
							array2D[i, j] = value;
						}
					}
				} else {
					int[] mIndex = new int[array.Rank];
					for (int index = 0; index < Count; index++) {
						StringUtil.LinearIndexToMultidimensionalIndex(index, strides, mIndex);
						array.SetValue(value, mIndex);
					}
				}
			} else {
				//ModifyAll(delegate(T item) { return (T)value.Clone(); });
				if (array1D != null) {
					for (int i = 0; i < array1D.Length; i++) {
						array1D[i] = (T)value.Clone();
					}
				} else if (array2D != null) {
					for (int i = 0; i < array2D.GetLength(0); i++) {
						for (int j = 0; j < array2D.GetLength(1); j++) {
							array2D[i, j] = (T)value.Clone();
						}
					}
				} else {
					int[] mIndex = new int[array.Rank];
					for (int index = 0; index < Count; index++) {
						StringUtil.LinearIndexToMultidimensionalIndex(index, strides, mIndex);
						array.SetValue((T)value.Clone(), mIndex);
					}
				}
			}
		}

		public void InitializeTo(DistributionArray<T> that)
		{
			if (IsValueType) {
				//ModifyAll(that.array, delegate(T item, T thatItem) { return thatItem; });
				if (array1D != null) {
					for (int i = 0; i < array1D.Length; i++) {
						array1D[i] = that.array1D[i];
					}
				} else if (array2D != null) {
					for (int i = 0; i < array2D.GetLength(0); i++) {
						for (int j = 0; j < array2D.GetLength(1); j++) {
							array2D[i, j] = that.array2D[i, j];
						}
					}
				} else {
					int[] mIndex = new int[array.Rank];
					for (int index = 0; index < Count; index++) {
						StringUtil.LinearIndexToMultidimensionalIndex(index, strides, mIndex);
						T item = (T)that.array.GetValue(mIndex);
						array.SetValue(item, mIndex);
					}
				}
			} else {
				//ModifyAll(that.array, delegate(T item, T thatItem) { return (T)thatItem.Clone(); });
				if (array1D != null) {
					for (int i = 0; i < array1D.Length; i++) {
						array1D[i] = (T)that.array1D[i].Clone();
					}
				} else if (array2D != null) {
					for (int i = 0; i < array2D.GetLength(0); i++) {
						for (int j = 0; j < array2D.GetLength(1); j++) {
							array2D[i, j] = (T)that.array2D[i, j].Clone();
						}
					}
				} else {
					int[] mIndex = new int[array.Rank];
					for (int index = 0; index < Count; index++) {
						StringUtil.LinearIndexToMultidimensionalIndex(index, strides, mIndex);
						T item = (T)that.array.GetValue(mIndex);
						array.SetValue((T)item.Clone(), mIndex);
					}
				}
			}
		}

		/// <summary>
		/// Set all distributions in the array to uniform
		/// </summary>
		public void SetToUniform()
		{
			ModifyAll(delegate(T item) { item.SetToUniform(); return item; });
		}

		/// <summary>
		/// Asks whether all distributions in a list are uniform
		/// </summary>
		/// <returns>True if all uniform. GFalse otherwise</returns>
		public bool IsUniform()
		{
			return Enumerable.TrueForAll<T>(this, delegate(T item) { return item.IsUniform(); });
		}

		/// <summary>
		/// Sets all distributions in the list to a specific distribution
		/// </summary>
		/// <param name="value">The distribution</param>
		public void SetTo(T value)
		{
			if (array1D != null) {
				for (int i = 0; i < array1D.Length; i++) {
					T item = array1D[i];
					item.SetTo(value);
					array1D[i] = item;
				}
			} else if (array2D != null) {
				for (int i = 0; i < array2D.GetLength(0); i++) {
					for (int j = 0; j < array2D.GetLength(1); j++) {
						T item = array2D[i, j];
						item.SetTo(value);
						array2D[i, j] = item;
					}
				}
			} else {
				int[] mIndex = new int[array.Rank];
				for (int index = 0; index < Count; index++) {
					StringUtil.LinearIndexToMultidimensionalIndex(index, strides, mIndex);
					T item = (T)array.GetValue(mIndex);
					item.SetTo(value);
					array.SetValue(item, mIndex);
				}
			}
		}

		/// <summary>
		/// Sets the current array of distributions to the given array of distributions
		/// </summary>
		/// <param name="value"></param>
		public void SetTo(DistributionArray<T> value)
		{
			if (IsValueType) {
				// Array.Copy does not work properly with reference types.
				if (value.Count != this.Count) throw new ArgumentException("value.Count (" + value.Count + ") != this.Count (" + this.Count + ")");
				Array.Copy(value.array, array, array.Length);
			} else {
				ModifyAll(value.array, delegate(T item, T aItem)
				{
					if (object.ReferenceEquals(item, null)) return (T)aItem.Clone();
					else item.SetTo(aItem); return item;
				});
			}
		}

		/// <summary>
		/// Clones the array and all items in the array.
		/// </summary>
		/// <returns>A new DistributionArray.</returns>
		public object Clone()
		{
			Array newArray = (Array)array.Clone();
			DistributionArray<T> result = new DistributionArray<T>(newArray);
			result.ModifyAll(delegate(T item) { return (T)item.Clone(); });
			return result;
		}

		/// <summary>
		/// Sets the current instance to an array of distributions each element
		/// of which is a product of the corresponding distributions in two given
		/// distribution arrays
		/// </summary>
		/// <param name="a">The first distribution array</param>
		/// <param name="b">The second distribution array</param>
		public void SetToProduct(DistributionArray<T> a, DistributionArray<T> b)
		{
			ModifyAll(a.array, b.array,
				delegate(T item, T aItem, T bItem)
				{
					item.SetToProduct(aItem, bItem);
					return item;
				});
		}

		/// <summary>
		/// Sets the current instance to an array of distributions each element
		/// of which is a ratio of the corrresponding distributions in two given
		/// distribution arrays
		/// </summary>
		/// <param name="numerator">The numerator distribution array</param>
		/// <param name="denominator">The denominator distribution array</param>
		public void SetToRatio(DistributionArray<T> numerator, DistributionArray<T> denominator)
		{
			ModifyAll(numerator.array, denominator.array,
				delegate(T item, T aItem, T bItem)
				{
					item.SetToRatio(aItem, bItem);
					return item;
				});
		}

		/// <summary>
		/// Sets the current instance to an array of distributions each element
		/// of which is a power of the corresponding element in a source distribution array
		/// </summary>
		/// <param name="a">The source distribution array</param>
		/// <param name="exponent">The exponent</param>
		public void SetToPower(DistributionArray<T> a, double exponent)
		{
			ModifyAll(a.array,
				delegate(T item, T aItem)
				{
					item.SetToPower(aItem, exponent);
					return item;
				});
		}

		/// <summary>
		/// Sets the current instance to an array of distributions each element
		/// of which is a weighted sum of the corresponding distributions in two given
		/// distribution arrays
		/// </summary>
		/// <param name="weight1">The first weight</param>
		/// <param name="a">The first distribution array</param>
		/// <param name="weight2">The second weight</param>
		/// <param name="b">The second distribution array</param>
		public void SetToSum(double weight1, DistributionArray<T> a, double weight2, DistributionArray<T> b)
		{
			ModifyAll(a.array, b.array,
				delegate(T item, T aItem, T bItem)
				{
					item.SetToSum(weight1, aItem, weight2, bItem);
					return item;
				});
		}

		/// <summary>
		/// The maximum difference across all corresponding distributions in
		/// this distribution array and that distribution array
		/// </summary>
		/// <param name="that">That distribution array</param>
		/// <returns>The maximum difference</returns>
		public double MaxDiff(object that)
		{
			DistributionArray<T> thatd = that as DistributionArray<T>;
			if ((object)thatd == null) return Double.PositiveInfinity;
			try {
				double diff = 0.0;
				ForEach(thatd.array,
					delegate(T item, T bItem)
					{
						diff = Math.Max(diff, item.MaxDiff(bItem));
					});
				return diff;
			} catch {
				return Double.PositiveInfinity;
			}
		}

		/// <summary>
		/// Ovverides the Equals method
		/// </summary>
		/// <param name="obj">The distribution array to compare to</param>
		/// <returns>true if the distribution arrays have equal values</returns>
		public override bool Equals(object obj)
		{
			return (MaxDiff(obj) == 0.0);
		}

		/// <summary>
		/// Overrides GetHashCode
		/// </summary>
		/// <returns>Hash code for the distribution array</returns>
		public override int GetHashCode()
		{
			int hash = Hash.Start;
			for (int dimension = 0; dimension < array.Rank; dimension++) {
				hash = Hash.Combine(hash, GetLength(dimension));
			}
			ForEach(delegate(T item) { hash = Hash.Combine(hash, item.GetHashCode()); });
			return hash;
		}

		/// <summary>
		/// The log-probability that two distributions would draw the same sample.
		/// </summary>
		/// <param name="that"></param>
		/// <returns><c>Math.Log(sum_x this.Evaluate(x)*that.Evaluate(x))</c></returns>
		/// <remarks>This can be considered a type of inner product between distributions.
		/// Another name might be "LogAverage" to go with "GetAverageLog".
		/// For a DistributionArray, this specializes to:
		/// <c>sum_i Math.Log(sum_x this[i].Evaluate(x)*that[i].Evaluate(x))</c>
		/// = <c>sum_i this[i].GetLogAverageOf(that[i])</c>
		/// </remarks>
		public double GetLogAverageOf(DistributionArray<T> that)
		{
			double sum = 0.0;
			ForEach(that.array, delegate(T item, T thatItem)
			{
				sum += item.GetLogAverageOf(thatItem);
			});
			return sum;
		}

		/// <summary>
		/// The expected logarithm of that distribution under this distribution.
		/// </summary>
		/// <param name="that">The distribution to take the logarithm of.</param>
		/// <returns><c>sum_x this.Evaluate(x)*Math.Log(that.Evaluate(x))</c></returns>
		/// <remarks>This is also known as the cross entropy.
		/// For a DistributionArray, this specializes to:
		/// <c>sum_i sum_x this[i].Evaluate(x)*Math.Log(that[i].Evaluate(x))</c>
		/// = <c>sum_i this[i].GetAverageLog(that[i])</c>
		/// </remarks>
		public double GetAverageLog(DistributionArray<T> that)
		{
			double sum = 0.0;
			ForEach(that.array, delegate(T item, T thatItem)
			{
				sum += item.GetAverageLog(thatItem);
			});
			return sum;
		}

		public static double EvaluateLn<DistributionType, DomainType>(DistributionArray<DistributionType> darray, DomainType[,] array)
		where DistributionType : CanGetLogProb<DomainType>,
			ICloneable,
			SettableTo<DistributionType>,
			SettableToUniform,
			Diffable,
			SettableToProduct<DistributionType>,
			SettableToRatio<DistributionType>,
			SettableToPower<DistributionType>,
			SettableToWeightedSum<DistributionType>,
			CanGetLogAverageOf<DistributionType>,
			CanGetAverageLog<DistributionType>
		{
			if (darray.array2D == null) throw new ArgumentException("DistributionArray has rank " + darray.array.Rank + ", point has rank " + array.Rank);
			if (darray.array2D.GetLength(0) != array.GetLength(0))
				throw new ArgumentException("DistributionArray lengths (" +
					StringUtil.CollectionToString(darray.GetLengths(), ",") + ") do not match array lengths (" +
					StringUtil.CollectionToString(StringUtil.ArrayDimensions(array), ",") + ")");
			double sum = 0.0;
			for (int i = 0; i < array.GetLength(0); i++) {
				for (int j = 0; j < array.GetLength(1); j++) {
					sum += darray[i, j].GetLogProb(array[i, j]);
				}
			}
			return sum;
		}
	}
#endif

#if false
	public class DistributionArray<T, DomainListType, DomainType> : DistributionArray<T>, IDistribution<DomainListType>
		//SettableToProduct<DistributionArray<T, DomainType>>,
		//SettableToRatio<DistributionArray<T, DomainType>>,
		//SettableTo<DistributionArray<T, DomainType>>
		//where DomainListType : IEnumerable<DomainType>
		where T : ICloneable, SettableTo<T>, SettableToUniform,
		Diffable,
		SettableToProduct<T>,
		SettableToRatio<T>,
		SettableToPower<T>,
		SettableToWeightedSum<T>,
		CanGetLogAverageOf<T>,
		CanGetAverageLog<T>,
		IDistribution<DomainType>
	{
		public DomainListType Point
		{
			get
			{
				throw new Exception("The method or operation is not implemented.");
			}
			set
			{
				throw new Exception("The method or operation is not implemented.");
			}
		}

		public bool IsPointMass
		{
			get { throw new Exception("The method or operation is not implemented."); }
		}

		public double GetLogProb(DomainListType value)
		{
			double sum = 0.0;
			if (value is Array) {
				ForEach<DomainType>((Array)(object)value, delegate(T distribution, DomainType item)
				{
					sum += distribution.GetLogProb(item);
				});
			} else {
				Enumerable.ForEach<T, DomainType>(this, (IEnumerable<DomainType>)value, delegate(T distribution, DomainType item)
				{
					sum += distribution.GetLogProb(item);
				});
			}
			return sum;
		}
	}
#endif

	internal interface ReducibleTo<T>
	{
		/// <summary>
		/// Remove dimensions via multiplication.
		/// </summary>
		/// <param name="keep">The dimensions to keep.</param>
		/// <param name="result">A distribution or distribution array.</param>
		/// <returns>An action which will perform the reduction.</returns>
		/// <remarks>Each element of result will be a product over the dimensions not kept.
		/// Result must already be the correct size.
		/// If keep is empty, no dimensions are kept so the result is a single distribution.
		/// Otherwise, result is a distribution array whose dimensions are the kept dimensions.</remarks>
		void ReduceTo(int[] keep, ICursorArray<T> result);
		void ReduceTo(T result);
	}

	/// <summary>
	/// The distribution of an array of independent variables, or equivalently
	/// an array of distributions.
	/// </summary>
	/// <typeparam name="DistributionType"></typeparam>
	/// <typeparam name="DomainType"></typeparam>
	/// <remarks>
	/// This class supports all of the IDistribution methods, as well as 
	/// being a CursorArray<typeparamref name="DistributionType"/>.
	/// To support plates, it implements a ReduceTo method which removes 
	/// dimensions via multiplication.
	/// </remarks>
	// TODO: named dimensions (name can be any object, such as Plate)
	internal class DistributionCursorArray<DistributionType, DomainType> :
		CursorArray<DistributionType>, ICursor, IDistribution<DomainType[]>,
		SettableTo<DistributionCursorArray<DistributionType, DomainType>>,
		SettableToProduct<DistributionCursorArray<DistributionType, DomainType>>,
		ReducibleTo<DistributionType>
		where DistributionType : IDistribution<DomainType>, ICursor, SettableTo<DistributionType>, SettableToProduct<DistributionType>
	{
		#region IDistribution methods
		/// <summary>
		/// True if all elements are constant.
		/// </summary>
		public bool IsPointMass
		{
			get
			{
				foreach (DistributionType dist in this) {
					if (!dist.IsPointMass) return false;
				}
				return true;
			}
		}

		public void SetToConstant()
		{
			foreach (DistributionType dist in this) {
				dist.Point = dist.Point;
			}
		}

		public DomainType[] Point
		{
			get
			{
				DomainType[] result = new DomainType[Count];
				int i = 0;
				foreach (DistributionType dist in this) {
					result[i++] = dist.Point;
				}
				// this only works if a distribution is stored sequentially in one array.
				// must be this[0] instead of cursor, since the Start property is important
				//CursorArray<ICursor> result = new CursorArray<ICursor>(this[0].Point, dim, stride);
				return result;
			}
			set
			{
				int i = 0;
				foreach (DistributionType dist in this) {
					dist.Point = value[i++];
				}
			}
		}

		public bool IsCompatibleWith(object thatd)
		{
			DistributionCursorArray<DistributionType, DomainType> that = thatd as DistributionCursorArray<DistributionType, DomainType>;
			if (that == null) return false;
			if (!base.IsCompatibleWith(that)) return false;
			return true;
			//return cursor.IsCompatibleWith(that.cursor);
		}
		public void CheckCompatible(object that)
		{
			if (!IsCompatibleWith(that))
				throw new ApplicationException("DistributionArrays are incompatible");
		}

		public void SetToUniform()
		{
			foreach (DistributionType d in this) {
				d.SetToUniform();
			}
		}
		public bool IsUniform()
		{
			foreach (DistributionType d in this) {
				if (!d.IsUniform()) return false;
			}
			return true;
		}

		public double GetLogProb(DomainType[] x)
		{
			throw new NotSupportedException();
		}

#if false
		public T[] Sample(T[] result)
		{
			Assert.IsTrue(result.Length == Count);
			int i = 0;
			foreach (DistributionType d in this) {
				result[i++] = d.Sample();
			}
			return result;
		}
		public T[] Sample()
		{
			return Sample(new T[Count]);
		}
#endif

		public void SetTo(DistributionCursorArray<DistributionType, DomainType> that)
		{
			CheckCompatible(that);
			Action action = delegate() { cursor.SetTo(that.cursor); };
			ForEach(that, action);
		}

		// a and b may be the same object as this.
		public void SetToProduct(DistributionCursorArray<DistributionType, DomainType> a, DistributionCursorArray<DistributionType, DomainType> b)
		{
			CheckCompatible(a);
			CheckCompatible(b);
			Action action = delegate() { cursor.SetToProduct(a.cursor, b.cursor); };
			ForEach(a, b, action);
		}

		public override bool Equals(object thatd)
		{
			DistributionCursorArray<DistributionType, DomainType> that = thatd as DistributionCursorArray<DistributionType, DomainType>;
			if (that == null) return false;
			if (!IsCompatibleWith(that)) return false;
			IEnumerator<DistributionType> iter = that.GetEnumerator();
			foreach (DistributionType d in this) {
				bool ok = iter.MoveNext();
				Assert.IsTrue(ok);
				if (!d.Equals(iter.Current)) return false;
			}
			return true;
		}
		public override int GetHashCode()
		{
			int hash = Hash.Start;
			for (int i = 0; i < Rank; i++) hash = Hash.Combine(hash, dim[i]);
			foreach (DistributionType d in this) {
				hash = Hash.Combine(hash, d.GetHashCode());
			}
			return hash;
		}
		public double MaxDiff(object thatd)
		{
			DistributionCursorArray<DistributionType, DomainType> that = thatd as DistributionCursorArray<DistributionType, DomainType>;
			if (that == null) return Double.PositiveInfinity;
			if (!IsCompatibleWith(that)) return Double.PositiveInfinity;
			double max = 0;
			IEnumerator<DistributionType> iter = that.GetEnumerator();
			foreach (DistributionType d in this) {
				bool ok = iter.MoveNext();
				Assert.IsTrue(ok);
				double diff = ((Diffable)d).MaxDiff(iter.Current);
				if (diff > max) max = diff;
			}
			return max;
		}
		#endregion

		/// <summary>
		/// Apply a reduction action to produce a smaller array.
		/// </summary>
		public void ReduceTo(int[] keep, ICursorArray<DistributionType> result)
		{
			Assert.IsTrue(keep.Length == result.Rank);
			int[] index = new int[Rank];
			int[] result_index = new int[result.Rank];
			for (int i = 0; i < Count; i++) {
				LinearIndexToMultidimensionalIndex(i, index);
				int index_sum = 0;
				for (int d = 0; d < index.Length; d++) {
					index_sum += index[d];
				}
				int result_index_sum = 0;
				for (int d = 0; d < keep.Length; d++) {
					result_index[d] = index[keep[d]];
					result_index_sum += result_index[d];
				}
				DistributionType d1 = this[index];
				DistributionType d2 = result[result_index];
				bool first = (result_index_sum == index_sum);
				if (first) {
					d2.SetTo(d1);
				} else {
					d2.SetToProduct(d2, d1);
				}
			}
		}
		/// <summary>
		/// Apply a reduction action to produce a single element distribution.
		/// </summary>
		public void ReduceTo(DistributionType result)
		{
			// FIXME: this could be more efficient
			bool first = true;
			foreach (DistributionType d in this) {
				if (first) {
					result.SetTo(d);
					first = false;
				} else {
					result.SetToProduct(result, d);
				}
			}
		}

		public DistributionCursorArray(DistributionType cursor, params int[] lengths)
			: base(cursor, lengths)
		{
		}

		public DistributionCursorArray(DistributionType cursor, IList<int> lengths, IList<int> stride)
			: base( cursor, lengths, stride)
		{
		}

		/// <summary>
		/// Add dimensions to an array by replication.
		/// </summary>
		/// <param name="lengths">The result array dimensions.</param>
		/// <param name="newPosition">For each original dimension d, newPosition[d] is its index in the 
		/// result dimensions.  Length == this.Rank.</param>
		/// <returns>A new array which uses the same storage but a different cursor.</returns>
		public DistributionCursorArray<DistributionType, DomainType> Replicate(int[] lengths, int[] newPosition)
		{
			int[] newStride = new int[lengths.Length]; // newStride[i] = 0
			for (int i = 0; i < Rank; i++) {
				newStride[newPosition[i]] = stride[i];
				Assert.IsTrue(lengths[newPosition[i]] == dim[i]);
			}
			return new DistributionCursorArray<DistributionType, DomainType>
				((DistributionType)this[0].ReferenceClone(), lengths, newStride);
		}

		/// <summary>
		/// Make a jagged array from a multidimensional array.
		/// </summary>
		/// <param name="isOuterDimension">For each original dimension d, 
		/// indicates whether it will be in the outer array. Length == this.Rank.</param>
		/// <returns>A jagged array [outer][inner] which uses the same storage but a different cursor.</returns>
		public DistributionCursorArray<DistributionCursorArray<DistributionType, DomainType>, DomainType[]>
						 Split(IList<bool> isOuterDimension)
		{
			int outerRank = StringUtil.Sum(isOuterDimension);
			if (outerRank == 0) throw new ArgumentException("No outer dimensions chosen");
			int innerRank = Rank - outerRank;
			if (innerRank == 0) throw new ArgumentException("No inner dimensions left");
			int[] innerLengths = new int[innerRank];
			int[] innerStride = new int[innerRank];
			int[] outerLengths = new int[outerRank];
			int[] outerStride = new int[outerRank];
			int innerIndex = 0;
			int outerIndex = 0;
			int i = 0;
			foreach (bool isOuter in isOuterDimension) {
				if (isOuter) {
					outerLengths[outerIndex] = dim[i];
					outerStride[outerIndex] = stride[i];
					outerIndex++;
				} else {
					innerLengths[innerIndex] = dim[i];
					innerStride[innerIndex] = stride[i];
					innerIndex++;
				}
				i++;
			}
			DistributionCursorArray<DistributionType, DomainType> inner =
				new DistributionCursorArray<DistributionType, DomainType>
				((DistributionType)this[0].ReferenceClone(), innerLengths, innerStride);
			return new DistributionCursorArray<DistributionCursorArray<DistributionType, DomainType>, DomainType[]>
				(inner, outerLengths, outerStride);
		}

		#region Copying
		public override ICursor ReferenceClone()
		{
			// must use this[0] not cursor
			return new DistributionCursorArray<DistributionType, DomainType>((DistributionType)this[0].ReferenceClone(), dim, stride);
		}
		public override object Clone()
		{
			return new DistributionCursorArray<DistributionType, DomainType>((DistributionType)this[0].ReferenceClone(), dim);
		}
		#endregion

		/// <summary>
		/// Overrides ToString method
		/// </summary>
		/// <returns>String representation of instance</returns>
		public override string ToString()
		{
			return StringUtil.ArrayToString(this, count, new int[Rank]);
#if false
			StringBuilder s = new StringBuilder();
			int index = 0;
			int[] indexList = new int[Rank];
			foreach (DistributionType d in this) {
				LinearIndexToMultidimensionalIndex(index++, indexList);
				string lhs = "["+Util.CollectionToString<int>(indexList,",")+"] ";
				s.Append(Util.JoinColumns("",lhs,d.ToString()));
				s.Append(Environment.NewLine);
			}
			return s.ToString();
#endif
		}
	}
}
