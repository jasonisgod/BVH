using System; 
using System.Collections.Generic;
using System.Diagnostics;
using System.Drawing;
using System.Threading;

namespace aabb
{
    internal static class _
    {
        public static void Resize<T>(this List<T> list, int size, T value = default(T)) where T : new()
        {
            list.Clear();
            for (int i = 0; i < size; i++)
            {
                list.Add(value != null? value: new T());
            }
        }

        public static List<T> List<T>(int size = 0, T value = default(T)) where T : new()
        {
            List<T> list = new List<T>();
            list.Resize(size, value);
            return list;
        }

        public static List<T> List<T>(IEnumerable<T> list)
        {
            return new List<T>(list);
        }

        public static int Count<K, V>(this Dictionary<K, V> dict, K key)
        {
            return Convert.ToInt32(dict.ContainsKey(key));
        }

        /*
        public static void Swap<T>(this List<T> list1, List<T> list2)
        {
            List<T> temp = _.List<T>(list1);
            list1.Clear();
            list1.AddRange(list2);
            list2.Clear();
            list2.AddRange(temp);
        }
        */

    }
    
}
