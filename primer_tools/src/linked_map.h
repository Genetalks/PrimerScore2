/*
 * linked_map
 *
 * An efficient map aggregate container that preserves the insertion
 * order of associated key/value pairs.
 *
 */

#ifndef _PT_LINKED_MAP_H_
#define _PT_LINKED_MAP_H_

#include <map>
#include <unordered_map>
#include <algorithm>
#include <iterator>
#include <list>

namespace pt {
   template <class K, class T, class M>
   class linked_map_impl {
   public:
      typedef linked_map_impl<K, T, M> linked_map_t;

      typedef typename M::size_type size_type;
      typedef std::pair<const K, T> value_type;
      typedef K key_type;
      typedef T mapped_type;
      typedef std::list<value_type> list_type;

      typedef typename list_type::iterator iterator;
      typedef typename list_type::const_iterator const_iterator;

      linked_map_impl ()
      { }

      linked_map_impl (list_type& list)
      {
         for (value_type& v : list) {
            insert (v);
         }
      }

      virtual ~linked_map_impl ()
      { }

      linked_map_impl<K, T, M>& operator= (
            const linked_map_impl<K, T, M>& map)
      {
         for (value_type v : map) {
            insert (v);
         }
         return *this;
      }

      bool empty () const
      {
         return iter_map.empty ();
      }

      size_type size () const
      {
         return iter_map.size ();
      }

      size_type max_size () const
      {
         return iter_map.max_size ();
      }

      size_type bucket_count () const
      {
         return iter_map.bucket_count ();
      }

      iterator begin ()
      {
         return value_list.begin ();
      }

      const_iterator begin () const
      {
         return value_list.begin ();
      }

      const_iterator cbegin () const
      {
         return value_list.cbegin ();
      }

      iterator end ()
      {
         return value_list.end ();
      }

      const_iterator end () const
      {
         return value_list.end ();
      }

      const_iterator cend () const
      {
         return value_list.cend ();
      }

      mapped_type& operator[] (const key_type& key)
      {
         if (! has_key (key)) {
            insert (value_type (key, T ()));
         }

         return at (key);
      }

      mapped_type& operator[] (key_type&& key)
      {
         if (! has_key (key)) {
            insert (value_type (key, T ()));
         }

         return at (key);
      }

      mapped_type& at (const key_type& key)
      {
         return (* iter_map.at (key)).second;
      }

      const mapped_type& at (const key_type& key) const
      {
         return (* iter_map.at (key)).second;
      }

      iterator find (const key_type& key)
      {
         auto iter = iter_map.find (key);

         if (iter != iter_map.end ()) {
            return iter->second;

         } else {
            return value_list.end ();
         }
      }

      const_iterator find (const key_type& key) const
      {
         auto iter = iter_map.find (key);
         if (iter != iter_map.end()){
            return iter->second;
         }
         else{
            return cend();
         }
      }

      bool has_key (const key_type& key) const
      {
         return iter_map.find (key) != iter_map.cend ();
      }

      iterator insert (const key_type& key, const mapped_type& value)
      {
         return insert (value_type (key, value));
      }

      iterator insert (const value_type& value)
      {
         return insert(end (), value);
      }

      iterator insert (const_iterator position, const key_type& key,
            const mapped_type& value)
      {
         return insert (position, value_type (key, value));
      }

      iterator insert (const_iterator position, const value_type& value)
      {
         auto iter = find (value.first);
         if (iter != end ()) {
            iter->second = value.second;
         }
         else{
             iter = value_list.insert (position, value);
             iter_map [value.first] = iter;
         }

         return iter;
      }

      iterator erase (const_iterator position)
      {
         iter_map.erase (position->first);
         return value_list.erase (position);
      }

      size_type erase (const key_type& key)
      {
         auto iter = find (key);

         if (iter != end ()) {
            erase (iter);
            return 1;

         } else {
            return 0;
         }
      }

      iterator erase (const_iterator first, const_iterator last)
      {
         iterator iter = first;

         while (iter != last && iter != cend ()) {
            iter = erase (iter);
         }

         return iter;
      }

      void clear ()
      {
         iter_map.clear ();
         value_list.clear ();
      }

      void swap (linked_map_impl<K, T, M>& map)
      {
         value_list.swap (map.value_list);
         iter_map.swap (map.iter_map);
      }

   private:
      M iter_map;
      list_type value_list;
   };

   template <class K, class T, class H = std::hash<K>, class C = std::equal_to<K>>
   using linked_hash = linked_map_impl<K, T,
      std::unordered_map<K, typename std::list<std::pair<const K, T>>::iterator, H, C> >;

   template <class K, class T, class C = std::less<K>>
   using linked_map = linked_map_impl<K, T,
      std::map<K,
                        typename std::list<std::pair<const K, T>>::iterator>>;
}

#endif //_GA_LINKED_MAP_H_
