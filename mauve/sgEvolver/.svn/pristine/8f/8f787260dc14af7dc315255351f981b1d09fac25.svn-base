#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#ifndef __GappedIntervalSequence_h__
#define __GappedIntervalSequence_h__

#include <vector>
#include "libGenome/IntervalSequenceTree.h"

/**
 * class implementing an Interval Sequence Tree
 * this is a tree for storing a changing sequence of intervals 
 * that was invented one rainy afternoon by Aaron Darling and
 * Michael Rusch.  
 * Important features are that stabbing queries,
 * stabbing insertions, and stabbing deletions are O( log n )
 * assuming a uniform distribution of stab sites.
 */
template< class Key, class Allocator = std::allocator<Key> >
class GappedIntervalSequence
{
public:
	typedef Key value_type;
//	typedef unsigned long long size_type;

	typedef Allocator                             allocator_type;
	typedef typename Allocator::reference         reference;
	typedef typename Allocator::const_reference   const_reference;
	typedef typename Allocator::size_type         size_type;
	typedef typename Allocator::difference_type   difference_type;
	typedef typename Allocator::pointer           pointer;
	typedef typename Allocator::const_pointer     const_pointer;

// node types and iterator definitions
protected:
	/**
	 * This class represents nodes of an Gapped Interval Sequence Tree.  Internal
	 * nodes define left and right to be non-null and key to
	 * be null.  Leaf nodes define left and right as null and key points
	 * to an interval.  The length field in an internal node is always the sum
	 * of lengths of the leaf nodes in its subtree.  The seq_length field is
	 * the actual length of sequence (not including gaps).
	 * The subtree_leaves field is
	 * defined as the number of leaf nodes below the node.
	 */
	class GisNode {
	public:			
		GisNode* parent;	/**< parent node pointer */
		GisNode* left;		/**< left node pointer */
		GisNode* right;		/**< right node pointer */
		size_type subtree_leaves;	/**< number of leaf nodes below this node */
		size_type length;		/**< total length of intervals below this node (sequence and gaps) */
		size_type seq_length;	/**< length of sequence below this node */
		Key* key;
		const GappedIntervalSequence< Key, Allocator >* gis_tree;
		GisNode( const GappedIntervalSequence< Key, Allocator >* gis_tree ) :
			parent( NULL ),
			left( NULL ),
			right( NULL ),
			subtree_leaves( 0 ),
			length( 0 ),
			seq_length( 0 ),
			key( NULL ),
			gis_tree( gis_tree )
			 {}
	};
	typedef typename Allocator::template rebind<GisNode>::other node_allocator_type;
	typedef typename node_allocator_type::pointer node_pointer;
	typedef typename node_allocator_type::const_pointer const_node_pointer;

//	typedef typename GisNode* node_pointer;
//	typedef typename const GisNode* const_node_pointer;

public:	

	// generic bidirectional iterator interface ripped from MSL, thanks guys
	template <bool is_const>
	class __generic_iterator
	{
	public:
		typedef typename GappedIntervalSequence::value_type       value_type;
//		typedef typename GappedIntervalSequence::difference_type  difference_type;
		typedef typename type_select<is_const, typename GappedIntervalSequence::const_pointer,
		                                  typename GappedIntervalSequence::pointer>::type pointer;
		typedef typename type_select<is_const, typename GappedIntervalSequence::const_reference,
		                                  typename GappedIntervalSequence::reference>::type reference;
		typedef bidirectional_iterator_tag        iterator_category;
		
		__generic_iterator() {}
		__generic_iterator(const __generic_iterator<false>& i) : ptr_(i.ptr_) {}
		reference operator * () const {return *ptr_->key;}
		pointer operator -> () const  {return ptr_->key;}
		__generic_iterator& operator ++ () {increment((const GisNode*&)ptr_); return *this;}
		__generic_iterator operator ++ (int) {__generic_iterator tmp(*this); ++(*this); return tmp;}
		__generic_iterator& operator -- () {decrement((const GisNode*&)ptr_); return *this;}
		__generic_iterator operator -- (int) {__generic_iterator tmp(*this); --(*this); return tmp;}
		friend bool operator ==(const __generic_iterator& x, const __generic_iterator& y) {return x.ptr_ == y.ptr_;}
		friend bool operator !=(const __generic_iterator& x, const __generic_iterator& y) {return x.ptr_ != y.ptr_;}
	private:
		typedef typename type_select<is_const, typename GappedIntervalSequence::const_node_pointer,
		                                  typename GappedIntervalSequence::node_pointer>::type node_pointer;

		node_pointer ptr_;

		explicit __generic_iterator(node_pointer n) : ptr_(n) {}

//		friend class __generic_iterator<true>;
		friend class GappedIntervalSequence;
	};

	friend class __generic_iterator<false>;
	friend class __generic_iterator<true>;
	typedef __generic_iterator<false> iterator;
	typedef __generic_iterator<true>  const_iterator;
	typedef std::reverse_iterator< iterator > reverse_iterator;
	typedef std::reverse_iterator< const_iterator > const_reverse_iterator;


// constructor related methods
	GappedIntervalSequence();
	template< class InputIterator >
	GappedIntervalSequence( InputIterator first, InputIterator last );
	GappedIntervalSequence( const GappedIntervalSequence& ist );
	GappedIntervalSequence& operator=( const GappedIntervalSequence& ist );
	~GappedIntervalSequence();

// standard container methods
	iterator begin();
	const_iterator begin() const;
	iterator end();
	const_iterator end() const;
	reverse_iterator rbegin();
	const_reverse_iterator rbegin() const;
	reverse_iterator rend();
	const_reverse_iterator rend() const;
	size_type max_size() const;
	bool empty() const;
	
// insertion and erasure
	iterator insert( const value_type& val, size_type point = IST_END );
	template <class InputIterator>
	void insert(InputIterator first, InputIterator last, size_type point = IST_END );

	size_type erase( size_type point, size_type length );
	void erase( iterator first, iterator last );

// search
	iterator find( size_type point );
	const_iterator find( size_type point ) const;

// interval sequence specific:
	/**
	 * Returns the total length of intervals contained in this interval sequence
	 */
	size_type length() const;
	size_type nodeCount() const;
	size_type countNodes( GisNode* x = NULL ) const;

	size_type sequenceLength() const;
	iterator find_seqindex( size_type seq_point );
	const_iterator find_seqindex( size_type seq_point ) const;
	
	size_type getSequenceStart( const_iterator iter ) const;
	size_type getStart( const_iterator iter ) const;
	size_type getSequenceStart( iterator iter ) const;
	size_type getStart( iterator iter ) const;
	
	void checkTree() const;
protected:
	GisNode *root;		/**< Root of the tree */
	GisNode *leftmost;	/**< Left most tree node, for begin() method */
	GisNode *rightmost;	/**< Right most tree node, for end() method */
	
	GisNode end_node;	/**< placeholder for the end node with a pointer to the root */
	
	static void propogateChanges( GisNode* cur_node, int64 length_diff, int64 seq_len_diff, int64 subtree_diff );
	const GisNode* recursiveFind( size_type& point, const GisNode* node ) const;
	const GisNode* recursiveSeqFind( size_type& point, const GisNode* const node ) const;
	static void increment( const GisNode*& x);
	static void increment( GisNode*& x);
	static void decrement( const GisNode*& x);
	static void decrement( GisNode*& x);
	static void deleteSubtree( GisNode*& istn );
	size_type recursiveGetSequenceStart( const GisNode* node ) const;
	size_type recursiveGetStart( const GisNode* node ) const;
	static void checkTree( GisNode* cur_node );
	static size_type checkNodeLengths( GisNode* node );
	static size_type checkNodeSeqLengths( GisNode* node );
	static bool isEndNode( const GisNode* node );
};



//template< class Key, class Allocator >
//GappedIntervalSequence< Key, Allocator >::IST_END = -1;

template< class Key, class Allocator >
inline
GappedIntervalSequence< Key, Allocator >::GappedIntervalSequence() :
root( NULL ),
leftmost( NULL ),
rightmost( NULL ),
end_node( this )
{
//	IST_END = -1;	// wraps to UINT64_MAX because IST_END is unsigned
}

template< class Key, class Allocator >
template< class InputIterator >
GappedIntervalSequence< Key, Allocator >::GappedIntervalSequence( InputIterator first, InputIterator last ) :
root( NULL ),
leftmost( NULL ),
rightmost( NULL ),
end_node( this )
{
	insert( first, last );
}

template< class Key, class Allocator >
GappedIntervalSequence< Key, Allocator >::GappedIntervalSequence( const GappedIntervalSequence& ist ) :
root( NULL ),
leftmost( NULL ),
rightmost( NULL ),
end_node( this )
{
	insert( ist.begin(), ist.end() );
}

template< class Key, class Allocator >
GappedIntervalSequence< Key, Allocator >& GappedIntervalSequence< Key, Allocator >::operator=( const GappedIntervalSequence& ist ){
	root = NULL;
	leftmost = NULL;
	rightmost = NULL;
	insert( ist.begin(), ist.end() );
	end_node = GisNode( this );
	return *this;
}

template< class Key, class Allocator >
GappedIntervalSequence< Key, Allocator >::~GappedIntervalSequence(){
	deleteSubtree( root );
}

template< class Key, class Allocator >
typename GappedIntervalSequence< Key, Allocator >::iterator 
GappedIntervalSequence< Key, Allocator >::begin(){
	GisNode* node = root;
	size_type point = 0;
	node = (GisNode*)recursiveFind( point, node );
	return iterator( node );
}

template< class Key, class Allocator >
typename GappedIntervalSequence< Key, Allocator >::const_iterator 
GappedIntervalSequence< Key, Allocator >::begin() const{
	const GisNode* node = root;
	size_type point = 0;
	node = recursiveFind( point, node );
	return const_iterator( node );
}

template< class Key, class Allocator >
typename GappedIntervalSequence< Key, Allocator >::iterator 
GappedIntervalSequence< Key, Allocator >::end(){
	return iterator( &end_node );
}

template< class Key, class Allocator >
typename GappedIntervalSequence< Key, Allocator >::const_iterator 
GappedIntervalSequence< Key, Allocator >::end() const{
	return const_iterator( &end_node );
}

template< class Key, class Allocator >
typename GappedIntervalSequence< Key, Allocator >::reverse_iterator 
GappedIntervalSequence< Key, Allocator >::rbegin(){
	return reverse_iterator( end() );
}

template< class Key, class Allocator >
typename GappedIntervalSequence< Key, Allocator >::const_reverse_iterator 
GappedIntervalSequence< Key, Allocator >::rbegin() const{
	return const_reverse_iterator( end() );
}

template< class Key, class Allocator >
typename GappedIntervalSequence< Key, Allocator >::reverse_iterator 
GappedIntervalSequence< Key, Allocator >::rend(){
	return reverse_iterator( begin() );
}

template< class Key, class Allocator >
typename GappedIntervalSequence< Key, Allocator >::const_reverse_iterator 
GappedIntervalSequence< Key, Allocator >::rend() const{
	return const_reverse_iterator( begin() );
}

template< class Key, class Allocator >
typename GappedIntervalSequence< Key, Allocator >::size_type 
GappedIntervalSequence< Key, Allocator >::max_size() const{
	return IST_END - 1;
}

template< class Key, class Allocator >
bool GappedIntervalSequence< Key, Allocator >::empty() const{
	return root == NULL ? true : false;
}

template< class Key, class Allocator >
void GappedIntervalSequence< Key, Allocator >::checkTree() const{
	checkTree( root );
	checkNodeLengths( root );
	checkNodeSeqLengths( root );
}

template< class Key, class Allocator >
void GappedIntervalSequence< Key, Allocator >::checkTree( 
	typename GappedIntervalSequence< Key, Allocator >::GisNode* cur_node ){
	if( cur_node ){
		if( cur_node->left && cur_node->left->parent != cur_node )
			cerr << "freakout\n";
		if( cur_node->right && cur_node->right->parent != cur_node )
			cerr << "freakout\n";
		checkTree( cur_node->left );
		checkTree( cur_node->right );
	}
}

template< class Key, class Allocator >
typename GappedIntervalSequence< Key, Allocator >::size_type
GappedIntervalSequence< Key, Allocator >::checkNodeLengths( 
	typename GappedIntervalSequence< Key, Allocator >::GisNode* cur_node ){
	if( cur_node ){
		if( cur_node->key ){
			if( cur_node->key->GetLength() != cur_node->length )
				cerr << "freakout\n";
			return cur_node->length;
		}
		size_type child_sum = checkNodeLengths( cur_node->left ) + checkNodeLengths( cur_node->right );
		if( child_sum != cur_node->length )
			cerr << "freakout\n";
		return cur_node->length;
	}
	return 0;
}

template< class Key, class Allocator >
typename GappedIntervalSequence< Key, Allocator >::size_type
GappedIntervalSequence< Key, Allocator >::checkNodeSeqLengths( 
	typename GappedIntervalSequence< Key, Allocator >::GisNode* cur_node ){
	if( cur_node ){
		if( cur_node->key ){
			if( cur_node->key->GetSeqLength() != cur_node->seq_length )
				cerr << "freakout\n";
			return cur_node->seq_length;
		}
		size_type child_sum = checkNodeSeqLengths( cur_node->left ) + checkNodeSeqLengths( cur_node->right );
		if( child_sum != cur_node->seq_length )
			cerr << "freakout\n";
		return cur_node->seq_length;
	}
	return 0;
}

template< class Key, class Allocator >
typename GappedIntervalSequence< Key, Allocator >::iterator 
GappedIntervalSequence< Key, Allocator >::insert( 
	const Key& val, 
	typename GappedIntervalSequence< Key, Allocator >::size_type point )
{
	size_type iv_offset = point;
	GisNode* ins_node = (GisNode*)recursiveFind( iv_offset, root );
	GisNode* new_node = new GisNode( this );
	new_node->key = new Key( val );
	new_node->length = val.GetLength();
	new_node->seq_length = val.GetSeqLength();
	new_node->subtree_leaves = 1;
	if( isEndNode( ins_node ) ){
		// end insert
		rightmost = new_node;
		if( root == NULL ){
			root = new_node;
			leftmost = new_node;
			return iterator( new_node );
		}
		// find the shallowest right insertion point
		ins_node = &end_node;
		decrement( ins_node );
		// make a new parent node
		GisNode* new_parent = new GisNode( this );
		new_parent->left = ins_node;
		new_parent->right = new_node;
		new_parent->parent = ins_node->parent;
		if( new_parent->parent == NULL )
			root = new_parent;
		else
			new_parent->parent->right = new_parent;
		ins_node->parent = new_parent;
		new_node->parent = new_parent;
		new_parent->length = ins_node->length;
		new_parent->seq_length = ins_node->seq_length;
		new_parent->subtree_leaves = 1;
		// update lengths and subtree_sizes along the path to the root
//		checkTree( root );
//		propogateChanges( new_node, 0, 0 );
//		propogateChanges( ins_node, 0, 0 );
		propogateChanges( new_parent, new_node->length, new_node->seq_length, 1 );
		return iterator( new_node );
	}

	// iv_offset is the distance into the node that the leaf should be split
	// 0 is a special case (left insert)
	if( iv_offset == 0 ){
		GisNode* new_parent = new GisNode( this );
		new_parent->left = new_node;
		new_parent->right = ins_node;
		new_parent->parent = ins_node->parent;
		if( new_parent->parent == NULL ){
			root = new_parent;
			leftmost = new_parent->left;
		}else if( new_parent->parent->right == ins_node )
			new_parent->parent->right = new_parent;
		else
			new_parent->parent->left = new_parent;
		new_parent->length = ins_node->length;
		new_parent->seq_length = ins_node->seq_length;
		new_parent->subtree_leaves = 1;

		ins_node->parent = new_parent;
		new_node->parent = new_parent;

		if( point == 0 )
			leftmost = new_node;
		// update lengths and subtree_sizes along the path to the root
//		checkTree( root );
//		propogateChanges( new_node, 0, 0 );
//		propogateChanges( ins_node, 0, 0 );
		propogateChanges( new_parent, new_node->length, new_node->seq_length, 1 );
	}else{
		// need to split a leaf node
		GisNode* new_gp = new GisNode( this );
		GisNode* new_parent = new GisNode( this );
		new_gp->parent = ins_node->parent;
		new_gp->right = new_parent;
		new_gp->left = new GisNode( this );
		new_gp->left->key = new Key( *ins_node->key );
		new_gp->left->key->CropEnd( ins_node->length - iv_offset );
		new_gp->left->length = new_gp->left->key->GetLength();
		new_gp->left->seq_length = new_gp->left->key->GetSeqLength();
		new_gp->left->subtree_leaves = 1;
		new_gp->left->parent = new_gp;
		
		ins_node->key->CropStart( iv_offset );
		ins_node->length = ins_node->key->GetLength();
		ins_node->seq_length = ins_node->key->GetSeqLength();
		ins_node->subtree_leaves = 1;
		ins_node->parent = new_parent;
		new_node->parent = new_parent;
		new_parent->left = new_node;
		new_parent->right = ins_node;
		new_parent->parent = new_gp;
		new_parent->length = new_node->length + ins_node->length;
		new_parent->seq_length = new_node->seq_length + ins_node->seq_length;
		new_parent->subtree_leaves = 2;

		new_gp->length = ins_node->length + new_gp->left->length;
		new_gp->seq_length = ins_node->seq_length + new_gp->left->seq_length;
		new_gp->subtree_leaves = 1;
		if( new_gp->parent == NULL ){
			root = new_gp;
			leftmost = new_gp->left;
			rightmost = ins_node;
		}else if( new_gp->parent->right == ins_node )
			new_gp->parent->right = new_gp;
		else
			new_gp->parent->left = new_gp;

		// update lengths and subtree_sizes along the path to the root
//		checkTree( root );
//		propogateChanges( new_node, 0, 0 );
//		propogateChanges( ins_node, 0, 0 );
//		propogateChanges( new_gp->left, 0, 0 );
//		propogateChanges( new_parent, 0, 0 );
		new_gp->subtree_leaves = 1;
		propogateChanges( new_gp, new_node->length, new_node->seq_length, 2 );
	}
	return iterator( new_node );
}

template< class Key, class Allocator >
template <class InputIterator>
void GappedIntervalSequence< Key, Allocator >::insert( 
	InputIterator first, 
	InputIterator last, 
	typename GappedIntervalSequence< Key, Allocator >::size_type point )
{

}

template< class Key, class Allocator >
typename GappedIntervalSequence< Key, Allocator >::size_type 
GappedIntervalSequence< Key, Allocator >::erase( 
	typename GappedIntervalSequence< Key, Allocator >::size_type point, 
	typename GappedIntervalSequence< Key, Allocator >::size_type length )
{
	size_type iv_offset = point;
	GisNode* ins_node = (GisNode*)recursiveFind( iv_offset, root );

	// iv_offset is the distance into the node that the leaf should be split
	// 0 is a special case (left delete)
	size_type deleted_nodes = 0;
	while( length > 0 ){
		if( isEndNode( ins_node ) ){
			// end delete?  that's illegal
			return deleted_nodes;
		}
		if( iv_offset == 0 ){
			if( length >= ins_node->length ){
				// delete the whole thing
				length -= ins_node->length;
				if( ins_node->parent == NULL ){
					// deleting the root
					delete ins_node;
					root = NULL;
					leftmost = NULL;
					rightmost = NULL;
					return deleted_nodes + 1;
				}

				GisNode* other_child = NULL, *del_node;
				if( ins_node->parent->left == ins_node ){
					other_child = ins_node->parent->right;
				}else if( ins_node->parent->right == ins_node ){
					other_child = ins_node->parent->left;
				}
				del_node = ins_node;
				increment( ins_node );

				// update tree structure
				GisNode* tmp_parent = other_child->parent;
				GisNode* tmp_gp = tmp_parent->parent;
				*tmp_parent = *other_child;
				tmp_parent->parent = tmp_gp;
				if( tmp_parent->left )
					tmp_parent->left->parent = tmp_parent;
				if( tmp_parent->right )
					tmp_parent->right->parent = tmp_parent;
				if( ins_node == other_child )
					ins_node = tmp_parent;
				delete other_child;
				
				// propogate deletion length thru root
//				checkTree( root );
				propogateChanges( tmp_gp, -del_node->length, -del_node->seq_length, -1 );
				if( del_node == leftmost )
					leftmost = tmp_parent;
				if( del_node == rightmost )
					rightmost = tmp_parent;
				// finally delete ins_node
				delete del_node;
				++deleted_nodes;
			}else{
				// crop from start
				int64 seq_len_diff = ins_node->key->GetSeqLength();
				ins_node->key->CropStart( length );
				seq_len_diff -= ins_node->key->GetSeqLength();
//				checkTree( root );
				propogateChanges( ins_node, -length, -seq_len_diff, 0 );
				return deleted_nodes;
			}
		}else if( length >= ins_node->length - iv_offset ){
			// crop from end
			int64 seq_len_diff = ins_node->key->GetSeqLength();
			ins_node->key->CropEnd( ins_node->length - iv_offset );
			seq_len_diff -= ins_node->key->GetSeqLength();
			length -= ins_node->length - iv_offset;
//			checkTree( root );
			propogateChanges( ins_node, -(ins_node->length - iv_offset), -seq_len_diff, 0 );
			increment( ins_node );
			iv_offset = 0;
		}else{
			// delete from middle (nastee part)
			int64 seq_len_diff = ins_node->seq_length;
			GisNode* new_parent = new GisNode( this );
			new_parent->left = ins_node;
			new_parent->length = ins_node->length;
			new_parent->seq_length = ins_node->seq_length;
			new_parent->subtree_leaves = 1;
			new_parent->right = new GisNode( this );
			new_parent->right->key = new Key( *ins_node->key );
			new_parent->right->length = ins_node->length - length - iv_offset;
			new_parent->right->key->CropStart( length + iv_offset );
			new_parent->right->seq_length = new_parent->right->key->GetSeqLength();
			new_parent->right->subtree_leaves = 1;
			new_parent->left->key->CropEnd( ins_node->length - iv_offset );
			new_parent->left->length = iv_offset;
			new_parent->left->seq_length = new_parent->left->key->GetSeqLength();
			new_parent->left->subtree_leaves = 1;
			new_parent->parent = ins_node->parent;
			seq_len_diff -= new_parent->left->seq_length;
			seq_len_diff -= new_parent->right->seq_length;
			if( new_parent->parent == NULL ){
				root = new_parent;
				rightmost = new_parent->right;
			}else if( new_parent->parent->left == ins_node )
				new_parent->parent->left = new_parent;
			else if( new_parent->parent->right == ins_node )
				new_parent->parent->right = new_parent;
			if( ins_node == rightmost )
				rightmost = new_parent->right;
			ins_node->parent = new_parent;
			new_parent->right->parent = new_parent;
//			checkTree( root );
//			propogateChanges( ins_node, 0, 0 );
//			propogateChanges( new_parent->right, 0, 0 );
			propogateChanges( new_parent, -length, -seq_len_diff, 1 );
			return deleted_nodes;
		}
	}
	return deleted_nodes;
}

template< class Key, class Allocator >
void GappedIntervalSequence< Key, Allocator >::propogateChanges( 
	GisNode* cur_node,
	int64 length_diff, 
	int64 seq_length_diff,
	int64 subtree_diff )
{
	vector< GisNode* > node_stack;
	while( cur_node != NULL ){
		if( cur_node->parent == cur_node )
			cerr << "when I say oh, you say shit!";
		cur_node->length += length_diff;
		cur_node->seq_length += seq_length_diff;
		cur_node->subtree_leaves += subtree_diff;
		node_stack.push_back( cur_node );
		cur_node = cur_node->parent;
	}
}

template< class Key, class Allocator >
void GappedIntervalSequence< Key, Allocator >::erase( 
	typename GappedIntervalSequence< Key, Allocator >::iterator first, 
	typename GappedIntervalSequence< Key, Allocator >::iterator last )
{

}

template< class Key, class Allocator >
typename GappedIntervalSequence< Key, Allocator >::iterator 
GappedIntervalSequence< Key, Allocator >::find( 
	typename GappedIntervalSequence< Key, Allocator >::size_type point ) 
{
	const GisNode* tmp = recursiveFind( point, root );
	return iterator( (GisNode*)tmp );
}

template< class Key, class Allocator >
typename GappedIntervalSequence< Key, Allocator >::const_iterator 
GappedIntervalSequence< Key, Allocator >::find( 
	typename GappedIntervalSequence< Key, Allocator >::size_type point ) const
{
	return const_iterator( GappedIntervalSequence< Key, Allocator >::recursiveFind( point, root ) );
}

template< class Key, class Allocator >
const typename GappedIntervalSequence< Key, Allocator >::GisNode* 
GappedIntervalSequence< Key, Allocator >::recursiveFind( 
	typename GappedIntervalSequence< Key, Allocator >::size_type& point, 
	const GisNode* node )  const{

	if( node == NULL )
		return &end_node;

	// return this node if it's a leaf
	if( node->key != NULL )
		return node;
	// look for the next node to recurse to
	if( point < node->length ){
		if( node->left ){
			if( point < node->left->length )
				return recursiveFind( point, node->left );
			point -= node->left->length;
		}
		return recursiveFind( point, node->right );
	}
	point -= node->length;
	// out of range
	return &(node->gis_tree->end_node);
}

template< class Key, class Allocator >
const typename GappedIntervalSequence< Key, Allocator >::GisNode* 
GappedIntervalSequence< Key, Allocator >::recursiveSeqFind( 
	typename GappedIntervalSequence< Key, Allocator >::size_type& seq_point, 
	const GisNode* const node )  const{

	if( node == NULL )
		return &end_node;

	// return this node if it's a leaf
	if( node->key != NULL )
		return node;
	// look for the next node to recurse to
	if( seq_point < node->length ){
		if( node->left ){
			if( seq_point < node->left->seq_length )
				return recursiveSeqFind( seq_point, node->left );
			seq_point -= node->left->seq_length;
		}
		return recursiveSeqFind( seq_point, node->right );
	}
	seq_point -= node->seq_length;
	// out of range
	return &(node->gis_tree->end_node);
}

template< class Key, class Allocator >
typename GappedIntervalSequence< Key, Allocator >::size_type 
GappedIntervalSequence< Key, Allocator >::sequenceLength() const {
	return root == NULL ? 0 : root->seq_length;
}

template< class Key, class Allocator >
typename GappedIntervalSequence< Key, Allocator >::iterator 
GappedIntervalSequence< Key, Allocator >::find_seqindex( 
	typename GappedIntervalSequence< Key, Allocator >::size_type seq_point ) 
{
	return iterator( recursiveSeqFind( seq_point, root ) );
}

template< class Key, class Allocator >
typename GappedIntervalSequence< Key, Allocator >::const_iterator 
GappedIntervalSequence< Key, Allocator >::find_seqindex( 
	typename GappedIntervalSequence< Key, Allocator >::size_type seq_point ) const {
	return const_iterator( recursiveSeqFind( seq_point, root ) );
}

template< class Key, class Allocator >
typename GappedIntervalSequence< Key, Allocator >::size_type
GappedIntervalSequence< Key, Allocator >::recursiveGetSequenceStart( const GisNode* node ) const{
	if( node == NULL || node == &end_node )
		return sequenceLength();
	if( node->parent == NULL )
		return 0;
	
	if( node->parent->right == node )
		return node->parent->left->seq_length + recursiveGetSequenceStart( node->parent );
	return recursiveGetSequenceStart( node->parent );
}

template< class Key, class Allocator >
typename GappedIntervalSequence< Key, Allocator >::size_type
GappedIntervalSequence< Key, Allocator >::recursiveGetStart( const GisNode* node )  const{
	if( node == NULL  || node == &end_node )
		return length();
	if( node->parent == NULL )
		return 0;
	
	if( node->parent->right == node )
		return node->parent->left->length + recursiveGetStart( node->parent );
	return recursiveGetStart( node->parent );
}

template< class Key, class Allocator >
typename GappedIntervalSequence< Key, Allocator >::size_type 
GappedIntervalSequence< Key, Allocator >::getSequenceStart( 
	typename GappedIntervalSequence< Key, Allocator >::const_iterator iter )  const
{
	return recursiveGetSequenceStart( iter.ptr_ );
}

template< class Key, class Allocator >
typename GappedIntervalSequence< Key, Allocator >::size_type 
GappedIntervalSequence< Key, Allocator >::getStart( 
	typename GappedIntervalSequence< Key, Allocator >::const_iterator iter )  const
{
	return recursiveGetStart( iter.ptr_ );
}

template< class Key, class Allocator >
typename GappedIntervalSequence< Key, Allocator >::size_type 
GappedIntervalSequence< Key, Allocator >::getSequenceStart( 
	typename GappedIntervalSequence< Key, Allocator >::iterator iter )  const
{
	return recursiveGetSequenceStart( iter.ptr_ );
}

template< class Key, class Allocator >
typename GappedIntervalSequence< Key, Allocator >::size_type 
GappedIntervalSequence< Key, Allocator >::getStart( 
	typename GappedIntervalSequence< Key, Allocator >::iterator iter )  const
{
	return recursiveGetStart( iter.ptr_ );
}

template< class Key, class Allocator >
typename GappedIntervalSequence< Key, Allocator >::size_type 
GappedIntervalSequence< Key, Allocator >::length() const{
	return root == NULL ? 0 : root->length;
}

template< class Key, class Allocator >
typename GappedIntervalSequence< Key, Allocator >::size_type 
GappedIntervalSequence< Key, Allocator >::nodeCount() const{
	return root == NULL ? 0 : root->subtree_leaves;
}

template< class Key, class Allocator >
typename GappedIntervalSequence< Key, Allocator >::size_type 
GappedIntervalSequence< Key, Allocator >::countNodes( GisNode* x ) const{
	if( x == NULL )
		x = root;
	if( x->key == NULL )
		return countNodes( x->left ) + countNodes( x->right ) + 1;
	return 1;
}


template< class Key, class Allocator >
void GappedIntervalSequence< Key, Allocator >::increment( const GisNode*& x){
	increment( (GisNode*&)x );
} 

template< class Key, class Allocator >
void GappedIntervalSequence< Key, Allocator >::decrement( const GisNode*& x) {
	decrement( (GisNode*&)x );
} 

template< class Key, class Allocator >
bool GappedIntervalSequence< Key, Allocator >::isEndNode( 
	const typename GappedIntervalSequence< Key, Allocator >::GisNode* x ) {
	if( x->parent == NULL && x->key == NULL && x->left == NULL )
		return true;
	return false;
}

template< class Key, class Allocator >
void GappedIntervalSequence< Key, Allocator >::increment( GisNode*& x){
	// find the least-ancestor with another child
	// and set x to that child
	while( x->parent != NULL ){
		if( x == x->parent->left &&
			x->parent->right != NULL ){
			x = x->parent->right;
			break;
		}else
			x = x->parent;
	}
	// if there were no other children to the right then we're at the end
	if( x->parent == NULL ){
		x = (GisNode*)&(x->gis_tree)->end_node;
		return;
	}

	// find the left-most leaf node below x
	while( x->key == NULL ){
		if( x->left != NULL )
			x = x->left;
		else if( x->right != NULL )
			x = x->right;
	}
}

template< class Key, class Allocator >
void GappedIntervalSequence< Key, Allocator >::decrement( GisNode*& x) {
	if( isEndNode( x ) ){
		x = x->gis_tree->root;	// x was the 'end' node
	}else{
		// find the least-ancestor with another child to the left
		// and set x to that child
		while( x->parent != NULL ){
			if( x == x->parent->right &&
				x->parent->left != NULL){
				x = x->parent->left;
				break;
			}else
				x = x->parent;
		}
		// if there was no other children to the left then we're at the end
		// raise hell! (cause an access violation)
		if( x->parent == NULL )
			x = NULL;
	}
	
	// find the right-most leaf node below x
	while( x->key == NULL ){
		if( x->right != NULL )
			x = x->right;
		else if( x->left != NULL )
			x = x->left;
	}
}

template< class Key, class Allocator >
void GappedIntervalSequence< Key, Allocator >::deleteSubtree( GisNode*& istn ) {
	if( istn == NULL )
		return;	// nothing more to delete
	deleteSubtree( istn->left );
	deleteSubtree( istn->right );
	if( istn->key != NULL )
		delete istn->key;
	delete istn;
}


#endif // __GappedIntervalSequence_h__
