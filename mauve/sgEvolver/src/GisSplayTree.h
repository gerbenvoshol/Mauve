#ifndef __GisSplayTree_h__
#define __GisSplayTree_h__

#include <vector>
#include "libGenome/IntervalSequenceTree.h"

/**
 * class implementing a Gapped Interval Sequence Tree
 * this tree stores a mapping between sequence and gapped sequence
 * coordinates in a potentially changing sequence
 * Important features are that stabbing queries,
 * stabbing insertions, and stabbing deletions are O( log n )
 */
template< class Key, class Allocator = std::allocator<Key> >
class GisSplayTree
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
	 */
	class GisNode {
	public:			
		GisNode* parent;	/**< parent node pointer */
		GisNode* left;		/**< left node pointer */
		GisNode* right;		/**< right node pointer */
		size_type length;		/**< total length of intervals below this node (sequence and gaps) */
		size_type seq_length;	/**< length of sequence below this node */
		Key* key;
		GisNode() :
			parent( NULL ),
			left( NULL ),
			right( NULL ),
			length( 0 ),
			seq_length( 0 ),
			key( NULL )
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
		typedef typename GisSplayTree::value_type       value_type;
//		typedef typename GisSplayTree::difference_type  difference_type;
		typedef typename type_select<is_const, typename GisSplayTree::const_pointer,
		                                  typename GisSplayTree::pointer>::type pointer;
		typedef typename type_select<is_const, typename GisSplayTree::const_reference,
		                                  typename GisSplayTree::reference>::type reference;
		typedef std::bidirectional_iterator_tag        iterator_category;
		
		__generic_iterator() {}
		__generic_iterator(const __generic_iterator<false>& i ) : ptr_(i.ptr_) {}
		reference operator * () const {return *ptr_->key;}
		pointer operator -> () const  {return ptr_->key;}
		__generic_iterator& operator ++ () {increment((const GisNode*&)ptr_); return *this;}
		__generic_iterator operator ++ (int) {__generic_iterator tmp(*this); ++(*this); return tmp;}
		__generic_iterator& operator -- () {decrement((const GisNode*&)ptr_); return *this;}
		__generic_iterator operator -- (int) {__generic_iterator tmp(*this); --(*this); return tmp;}
		friend bool operator ==(const __generic_iterator& x, const __generic_iterator& y) {return x.ptr_ == y.ptr_;}
		friend bool operator !=(const __generic_iterator& x, const __generic_iterator& y) {return x.ptr_ != y.ptr_;}
	private:
		typedef typename type_select<is_const, typename GisSplayTree::const_node_pointer,
		                                  typename GisSplayTree::node_pointer>::type node_pointer;

		node_pointer ptr_;

		explicit __generic_iterator(node_pointer n) : ptr_(n) {}

//		friend class __generic_iterator<true>;
		friend class GisSplayTree;
	};

	friend class __generic_iterator<false>;
	friend class __generic_iterator<true>;
	typedef __generic_iterator<false> iterator;
	typedef __generic_iterator<true>  const_iterator;
	typedef std::reverse_iterator< iterator > reverse_iterator;
	typedef std::reverse_iterator< const_iterator > const_reverse_iterator;


// constructor related methods
	GisSplayTree();
	template< class InputIterator >
	GisSplayTree( InputIterator first, InputIterator last );
	GisSplayTree( const GisSplayTree& ist );
	GisSplayTree& operator=( const GisSplayTree& ist );
	~GisSplayTree();

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
	
	size_type getSequenceStart( const_iterator iter );
	size_type getStart( const_iterator iter );
	size_type getSequenceStart( iterator iter );
	size_type getStart( iterator iter );
	
	void checkTree() const;

	static int debug_checking;
protected:
	GisNode *root;		/**< Root of the tree */
	
	GisNode* recursiveFind( size_type& point, GisNode* node ) const;
	GisNode* recursiveSeqFind( size_type& point, GisNode* node ) const;
	static void increment( const GisNode*& x);
	static void increment( GisNode*& x);
	void decrement( const GisNode*& x);
	void decrement( GisNode*& x);
	void deleteTree();
	size_type recursiveGetSequenceStart( const GisNode* node ) const;
	size_type recursiveGetStart( const GisNode* node ) const;
	static void checkTree( GisNode* cur_node );
	static size_type checkNodeLengths( GisNode* node );
	static size_type checkNodeSeqLengths( GisNode* node );

	void splay( GisNode* node );
	void recalculateLengths( GisNode* n );
};


template< class Key, class Allocator >
int GisSplayTree< Key, Allocator >::debug_checking = 0;


//template< class Key, class Allocator >
//GisSplayTree< Key, Allocator >::IST_END = -1;

template< class Key, class Allocator >
inline
GisSplayTree< Key, Allocator >::GisSplayTree() :
root( NULL )
{
}

template< class Key, class Allocator >
template< class InputIterator >
GisSplayTree< Key, Allocator >::GisSplayTree( InputIterator first, InputIterator last ) :
root( NULL )
{
	insert( first, last );
}

template< class Key, class Allocator >
GisSplayTree< Key, Allocator >::GisSplayTree( const GisSplayTree& ist ) :
root( NULL )
{
	insert( ist.begin(), ist.end() );
}

template< class Key, class Allocator >
GisSplayTree< Key, Allocator >& GisSplayTree< Key, Allocator >::operator=( const GisSplayTree& ist ){
	root = NULL;
	insert( ist.begin(), ist.end() );
	return *this;
}

template< class Key, class Allocator >
GisSplayTree< Key, Allocator >::~GisSplayTree(){
	deleteTree();
}

template< class Key, class Allocator >
typename GisSplayTree< Key, Allocator >::iterator 
GisSplayTree< Key, Allocator >::begin(){
	GisNode* node = root;
	size_type point = 0;
	node = (GisNode*)recursiveFind( point, node );
	return iterator( node );
}

template< class Key, class Allocator >
typename GisSplayTree< Key, Allocator >::const_iterator 
GisSplayTree< Key, Allocator >::begin() const{
	GisNode* node = root;
	size_type point = 0;
	node = recursiveFind( point, node );
	return const_iterator( node );
}

template< class Key, class Allocator >
typename GisSplayTree< Key, Allocator >::iterator 
GisSplayTree< Key, Allocator >::end(){
	GisNode* nowhere = NULL;
	return iterator( nowhere );
}

template< class Key, class Allocator >
typename GisSplayTree< Key, Allocator >::const_iterator 
GisSplayTree< Key, Allocator >::end() const{
	GisNode* nowhere = NULL;
	return const_iterator( nowhere );
}

template< class Key, class Allocator >
typename GisSplayTree< Key, Allocator >::reverse_iterator 
GisSplayTree< Key, Allocator >::rbegin(){
	return reverse_iterator( end() );
}

template< class Key, class Allocator >
typename GisSplayTree< Key, Allocator >::const_reverse_iterator 
GisSplayTree< Key, Allocator >::rbegin() const{
	return const_reverse_iterator( end() );
}

template< class Key, class Allocator >
typename GisSplayTree< Key, Allocator >::reverse_iterator 
GisSplayTree< Key, Allocator >::rend(){
	return reverse_iterator( begin() );
}

template< class Key, class Allocator >
typename GisSplayTree< Key, Allocator >::const_reverse_iterator 
GisSplayTree< Key, Allocator >::rend() const{
	return const_reverse_iterator( begin() );
}

template< class Key, class Allocator >
typename GisSplayTree< Key, Allocator >::size_type 
GisSplayTree< Key, Allocator >::max_size() const{
	return IST_END - 1;
}

template< class Key, class Allocator >
bool GisSplayTree< Key, Allocator >::empty() const{
	return root == NULL ? true : false;
}

template< class Key, class Allocator >
void GisSplayTree< Key, Allocator >::checkTree() const{
	if( debug_checking > 0 ){
		checkTree( root );
		checkNodeLengths( root );
		checkNodeSeqLengths( root );
	}
}

template< class Key, class Allocator >
void GisSplayTree< Key, Allocator >::checkTree( 
	typename GisSplayTree< Key, Allocator >::GisNode* cur_node ){
	if( cur_node ){
		if( cur_node->left && cur_node->left->parent != cur_node )
			std::cerr << "freakout\n";
		if( cur_node->right && cur_node->right->parent != cur_node )
			std::cerr << "freakout\n";
		checkTree( cur_node->left );
		checkTree( cur_node->right );
	}
}

/**
 * validate that recorded subtree lengths and node lengths are consistent
 */
template< class Key, class Allocator >
typename GisSplayTree< Key, Allocator >::size_type
GisSplayTree< Key, Allocator >::checkNodeLengths( 
	typename GisSplayTree< Key, Allocator >::GisNode* cur_node )
{
	if( cur_node == NULL )
		return 0;

	long left_len = cur_node->left != NULL ? cur_node->left->length : 0;
	long right_len = cur_node->right != NULL ? cur_node->right->length : 0;
	// do left_len and right_len match the actual subtree lengths?
	if( left_len != checkNodeLengths( cur_node->left ) )
		std::cerr << "freakout at:" << __FILE__ << ":" << __LINE__ << "\n";
	if( right_len != checkNodeLengths( cur_node->right ) )
		std::cerr << "freakout at:" << __FILE__ << ":" << __LINE__ << "\n";
	// do they all sum up to the correct value?
	if( left_len + right_len + cur_node->key->GetLength() != cur_node->length ){
		std::cerr << "freakout at:" << __FILE__ << ":" << __LINE__ << "\n";
		std::cerr << cur_node->key;
	}
	return cur_node->length;
}

	/**
	 * validate that the recorded sequence length and actual sequence length
	 * below a node are consistent with each other
	 */
template< class Key, class Allocator >
typename GisSplayTree< Key, Allocator >::size_type
GisSplayTree< Key, Allocator >::checkNodeSeqLengths( 
	typename GisSplayTree< Key, Allocator >::GisNode* cur_node )
{
	if( cur_node == NULL )
		return 0;

	long left_len = cur_node->left != NULL ? cur_node->left->seq_length : 0;
	long right_len = cur_node->right != NULL ? cur_node->right->seq_length : 0;
	// do left_len and right_len match the actual subtree lengths?
	if( left_len != checkNodeSeqLengths( cur_node->left ) )
		std::cerr << "freakout at:" << __FILE__ << ":" << __LINE__ << "\n";
	if( right_len != checkNodeSeqLengths( cur_node->right ) )
		std::cerr << "freakout at:" << __FILE__ << ":" << __LINE__ << "\n";
	// do they all sum up to the correct value?
	if( left_len + right_len + cur_node->key->GetSeqLength() != cur_node->seq_length ){
		std::cerr << "freakout at:" << __FILE__ << ":" << __LINE__ << "\n";
		std::cerr << cur_node->key;
	}
	return cur_node->seq_length;
}


template< class Key, class Allocator >
typename GisSplayTree< Key, Allocator >::iterator 
GisSplayTree< Key, Allocator >::insert( 
	const Key& val, 
	typename GisSplayTree< Key, Allocator >::size_type point )
{
	if( point == IST_END )
		point = root == NULL ? 0 : root->length;
	if( ( root == NULL && point > 0 ) || ( root != NULL && point > root->length ) ){
		if( root == NULL ){
			std::cerr << "root is null and point is: " << point << std::endl;
			point = 0;
		}else{
			point = root->length;
			std::cerr << "root->length is: " << root->length << ", point is: " << point << std::endl;
		}
		std::cerr << "Working around the problem\n";
//		throw "GisSplayTree::insert() out of bounds";
	}
	size_type position = point;
	GisNode* new_node = new GisNode();
	new_node->key = new Key( val );
	new_node->length = val.GetLength();
	new_node->seq_length = val.GetSeqLength();

	// just insert new_node as root if the tree is empty
	if( root == NULL ){
		root = new_node;
		//checkTree();
		return iterator( new_node );
	}

	GisNode* ins_node = recursiveFind( position, root );
	checkTree();

	// insert the new node below ins_node
	if( position > 0 && position < ins_node->key->GetLength() ){
		// trunc ins_node, do a right insert of new_node and the right part of ins_node
		GisNode* ins_right = new GisNode();
		// Question:  does inserting two nodes at once violate the splay rules?
		// probably, but i'm not sure it really matters
		new_node->right = ins_right;
		ins_right->parent = new_node;
		ins_right->key = new Key( *(ins_node->key) );
		// crop the key
		ins_right->key->CropStart( position );
		ins_right->length = ins_right->key->GetLength();
		ins_right->seq_length = ins_right->key->GetSeqLength();
		// question: do I need to update new_node.length or will the splay operations 
		// take care of it?

		ins_node->key->CropEnd( ins_node->key->GetLength() - position );
		ins_node->seq_length = ins_node->key->GetSeqLength();
		ins_node->length = ins_node->key->GetLength();
		// now position[0] ought to be equal to ins_node.length,
		// so new_node should get inserted right below it
	}

	if( position == 0 ){
		// find the right-most child of the left subtree and insert there
		GisNode* cur_node = ins_node->left;
		if( cur_node != NULL ){
			while( cur_node->right != NULL )
				cur_node = cur_node->right;
			// insert to the right of cur_node
			cur_node->right = new_node;
			new_node->parent = cur_node;
		}else{
			// insert to the left of ins_node
			ins_node->left = new_node;
			new_node->parent = ins_node;
		}
	}
	if( position >= ins_node->key->GetLength() ){
		// find the left-most child of the right subtree and insert there
		GisNode* cur_node = ins_node->right;
		if( cur_node != NULL ){
			while( cur_node->left != NULL )
				cur_node = cur_node->left;
			// insert to the left of cur_node
			cur_node->left = new_node;
			new_node->parent = cur_node;
		}else{
			// insert to the right of ins_node
			ins_node->right = new_node;
			new_node->parent = ins_node;
		}
	}
	splay( new_node );
	checkTree();
	return iterator( new_node );
}

template< class Key, class Allocator >
template <class InputIterator>
void GisSplayTree< Key, Allocator >::insert( 
	InputIterator first, 
	InputIterator last, 
	typename GisSplayTree< Key, Allocator >::size_type point )
{
	size_type cur_point = point;
	while( first != last ){
		insert( *(first.ptr_->key), point );
		point += first.ptr_->key->GetLength();
		first++;
	}
}

/**
 * remove a range of coordinates from the sequence
 */
template< class Key, class Allocator >
typename GisSplayTree< Key, Allocator >::size_type 
GisSplayTree< Key, Allocator >::erase( 
	typename GisSplayTree< Key, Allocator >::size_type point, 
	typename GisSplayTree< Key, Allocator >::size_type length )
{
	// loop while there are nodes to delete
	while( length > 0 ){
		// splay delete:  splay node X, splay node X-1
		//                delete X, patch X's right subtree onto X-1
		size_type position = point;
		GisNode* del_node = recursiveFind( position, root );
		splay( del_node );
		checkTree();

		GisNode* left_ins = NULL;
		if( position > 0 ){
			// just crop this node down, don't delete it
			long left_len = del_node->left != NULL ? del_node->left->length : 0;
			if( position + length < left_len + del_node->key->GetLength() ){
				// this node needs something deleted from its middle->->->
				// split it into two separate nodes at the left delete
				// point, insert the preserved left-side portion of the node
				// and then do a single delete operation next time thru the loop
				Key left_k( *(del_node->key) );
				left_k.CropEnd( left_k.GetLength() - position );
				del_node->key->CropStart( position );
				recalculateLengths( del_node );
				checkTree();
				insert( left_k, point - position );
			}else{
				// crop from right
				length -= del_node->key->GetLength() - position;
				del_node->key->CropEnd( del_node->key->GetLength() - position );
				recalculateLengths( del_node );
			}
		}else if( length < del_node->key->GetLength() ){
			// crop from the left, we're done
			del_node->key->CropStart( length );
			recalculateLengths( del_node );
			length = 0;
		}else{
			// delete the whole darn thing
			if( root->left != NULL ){
				position = root->left->length - 1;
				GisNode* prev_node = recursiveFind( position, root );
				splay( prev_node );
				checkTree();
				if( root->right->right != NULL )
					root->right->right->parent = root;
				root->right = root->right->right;
				root->length -= del_node->key->GetLength();
				root->seq_length -= del_node->key->GetSeqLength();
			}else{
				// just delete root
				if( root->right != NULL )
					root->right->parent = NULL;
				root = root->right;
			}

//			point += del_node->key->GetLength();
			length -= del_node->key->GetLength();
			delete del_node->key;
			delete del_node;
		}

	}
	return 0;
}

template< class Key, class Allocator >
void GisSplayTree< Key, Allocator >::erase( 
	typename GisSplayTree< Key, Allocator >::iterator first, 
	typename GisSplayTree< Key, Allocator >::iterator last )
{
	throw "Implement me!";
}

/**
 * splay a node to the root of the tree
 */
template< class Key, class Allocator >
void GisSplayTree< Key, Allocator >::splay( 
	typename GisSplayTree< Key, Allocator >::GisNode* node )
{ 
	// splay operations and node naming convention taken from
	// http://www.cs.nyu.edu/algvis/java/SplayTree.html

	while( node->parent != NULL ){
		GisNode* x = node;
		GisNode* y = x->parent;
		GisNode* z = y->parent != NULL ? y->parent : y;
		if( node->parent->left == node ){
			// zigging
			if( node->parent->parent == NULL ){
				// zig
				y->left = x->right;
				x->right = y;
			}else if( node->parent->parent->left == node->parent ){
				// zig-zig
				z->left = y->right;
				y->left = x->right;
				y->right = z;
				x->right = y;
			}else{
				// zag-zig
				z->right = x->left;
				y->left = x->right;
				x->right = y;
				x->left = z;
			}
		}else{
			if( node->parent->right != node )
				throw "Bad error on line 544";
			// zagging
			if( node->parent->parent == NULL ){
				// zag
				y->right = x->left;
				x->left = y;
			}else if( node->parent->parent->left == node->parent ){
				// zig-zag
				z->left = x->right;
				y->right = x->left;
				x->left = y;
				x->right = z;
			}else{
				// zag-zag
				z->right = y->left;
				y->right = x->left;
				y->left = z;
				x->left = y;
			}
		}
		// update parents and lengths
		x->parent = z->parent;
		if( x->parent != NULL ){
			if( x->parent->left == z )
				x->parent->left = x;
			else
				x->parent->right = x;
		}
		recalculateLengths( z );
		recalculateLengths( y );
		recalculateLengths( x );
	}
	root = node;
	//checkTree();
}



/**
 * recalculate a node's lengths and update its children's
 * parent pointers
 */
template< class Key, class Allocator >
void GisSplayTree< Key, Allocator >::recalculateLengths( 
	typename GisSplayTree< Key, Allocator >::GisNode* n )
{ 
	n->length = n->key->GetLength();
	n->length += n->right != NULL ? n->right->length : 0;
	n->length += n->left != NULL ? n->left->length : 0;
	n->seq_length = n->key->GetSeqLength();
	n->seq_length += n->right != NULL ? n->right->seq_length : 0;
	n->seq_length += n->left != NULL ? n->left->seq_length : 0;

	if( n->left != NULL )
		n->left->parent = n;
	if( n->right != NULL )
		n->right->parent = n;
}


/**
 * find the node below cur_node containing a given position in the gapped sequence,
 * starting at the left-most position below cur_node
 * @param position	The position to search for.  A single element array.  The value
 *                  is modified to reflect the distance into the returned node where
 *                  the requested position actually occurs.
 * @param cur_node	The tree node to use as the base for the search.  usually the root
 * version 2 of recursiveFind -- don't build a potentially huge call stack!
 */
template< class Key, class Allocator >
typename GisSplayTree< Key, Allocator >::GisNode* 
GisSplayTree< Key, Allocator >::recursiveFind( 
	typename GisSplayTree< Key, Allocator >::size_type& position, 
	GisNode* cur_node )  const{

	while( cur_node != NULL ){

		size_type left_len = cur_node->left != NULL ? cur_node->left->length : 0;
		if( cur_node->left != NULL && position < cur_node->left->length ){
			cur_node = cur_node->left;
			continue;
		}

		// it's not part of the left subtree, subtract off the left subtree length
		position -= left_len;
		if( cur_node->right != NULL && position >= cur_node->key->GetLength() ){
			position -= cur_node->key->GetLength();
			cur_node = cur_node->right;
			continue;
		}

		// return this node if nothing else can be done
		return cur_node;
	}
	return NULL;
}


/**
 * Find the node below cur_node containing a given position in the ungapped sequence,
 * starting at the left-most position below cur_node
 * @param position	The position to search for.  A single element array.  The value
 *                  is modified to reflect the distance into the returned node where
 *                  the requested position actually occurs.
 * @param cur_node	The tree node to use as the base for the search.  usually the root
 * version 2 of recursiveSeqFind -- don't build a potentially huge call stack!
 */
template< class Key, class Allocator >
typename GisSplayTree< Key, Allocator >::GisNode* 
GisSplayTree< Key, Allocator >::recursiveSeqFind( 
	typename GisSplayTree< Key, Allocator >::size_type& position, 
	GisNode* cur_node )  const{

	while( cur_node != NULL ){
		long left_len = cur_node->left != NULL ? cur_node->left->seq_length : 0;
		if( cur_node->left != NULL && position < cur_node->left->seq_length ){
			cur_node = cur_node->left;
			continue;
		}

		// it's not part of the left subtree, subtract off the left subtree length
		position -= left_len;
		if( cur_node->right != NULL && position >= cur_node->key->GetSeqLength() ){
			position -= cur_node->key->GetSeqLength();
			cur_node = cur_node->right;
			continue;
		}

		// return this node if nothing else can be done
		return cur_node;
	}
	return NULL;
}



/** find the interval containing a position in the gapped sequence */
template< class Key, class Allocator >
typename GisSplayTree< Key, Allocator >::iterator 
GisSplayTree< Key, Allocator >::find( 
	typename GisSplayTree< Key, Allocator >::size_type point ) 
{
	GisNode* node = recursiveFind( point, root );
	splay( node );
	return iterator( root );
}



template< class Key, class Allocator >
typename GisSplayTree< Key, Allocator >::const_iterator 
GisSplayTree< Key, Allocator >::find( 
	typename GisSplayTree< Key, Allocator >::size_type point ) const
{
	return const_iterator( GisSplayTree< Key, Allocator >::recursiveFind( point, root ) );
}



/** returns the length of ungapped sequence stored in the tree */
template< class Key, class Allocator >
typename GisSplayTree< Key, Allocator >::size_type 
GisSplayTree< Key, Allocator >::sequenceLength() const {
	return root == NULL ? 0 : root->seq_length;
}

/** find the interval containing a position in the ungapped sequence */
template< class Key, class Allocator >
typename GisSplayTree< Key, Allocator >::iterator 
GisSplayTree< Key, Allocator >::find_seqindex( 
	typename GisSplayTree< Key, Allocator >::size_type seq_point ) 
{
	GisNode* node = recursiveSeqFind( seq_point, root );
	splay( node );
	return iterator( root );
}

template< class Key, class Allocator >
typename GisSplayTree< Key, Allocator >::const_iterator 
GisSplayTree< Key, Allocator >::find_seqindex( 
	typename GisSplayTree< Key, Allocator >::size_type seq_point ) const
{
	return const_iterator( recursiveSeqFind( seq_point, root ) );
}

template< class Key, class Allocator >
typename GisSplayTree< Key, Allocator >::size_type 
GisSplayTree< Key, Allocator >::getSequenceStart( 
	typename GisSplayTree< Key, Allocator >::const_iterator iter )
{
	splay( iter.ptr_ );
	return root->left != NULL ? root->left->seq_length : 0;
}

template< class Key, class Allocator >
typename GisSplayTree< Key, Allocator >::size_type 
GisSplayTree< Key, Allocator >::getStart( 
	typename GisSplayTree< Key, Allocator >::const_iterator iter )
{
	splay( iter.ptr_ );
	return root->left != NULL ? root->left->length : 0;
}

template< class Key, class Allocator >
typename GisSplayTree< Key, Allocator >::size_type 
GisSplayTree< Key, Allocator >::getSequenceStart( 
	typename GisSplayTree< Key, Allocator >::iterator iter )
{
	splay( iter.ptr_ );
	return root->left != NULL ? root->left->seq_length : 0;
}

template< class Key, class Allocator >
typename GisSplayTree< Key, Allocator >::size_type 
GisSplayTree< Key, Allocator >::getStart( 
	typename GisSplayTree< Key, Allocator >::iterator iter )
{
	splay( iter.ptr_ );
	return root->left != NULL ? root->left->length : 0;
}

/** returns the total length of the gapped sequence stored in the tree */
template< class Key, class Allocator >
typename GisSplayTree< Key, Allocator >::size_type 
GisSplayTree< Key, Allocator >::length() const{
	return root == NULL ? 0 : root->length;
}

template< class Key, class Allocator >
typename GisSplayTree< Key, Allocator >::size_type 
GisSplayTree< Key, Allocator >::nodeCount() const{
	return root == NULL ? 0 : countNodes(root);
}

template< class Key, class Allocator >
typename GisSplayTree< Key, Allocator >::size_type 
GisSplayTree< Key, Allocator >::countNodes( GisNode* x ) const{
	if( x == NULL )
		return 0;
	return countNodes( x->left ) + countNodes( x->right ) + 1;
}


template< class Key, class Allocator >
void GisSplayTree< Key, Allocator >::increment( const GisNode*& x){
	increment( (GisNode*&)x );
} 

template< class Key, class Allocator >
void GisSplayTree< Key, Allocator >::decrement( const GisNode*& x) {
	decrement( (GisNode*&)x );
} 

template< class Key, class Allocator >
void GisSplayTree< Key, Allocator >::increment( GisNode*& x){
	if( x == NULL )
		// x is already the right-most tree value
		// throw an exception?
		throw "incrementing out of bounds!";

	// if x has a right child, find its leftmost descendant
	if( x->right != NULL ){
		x = x->right;
		while( x->left != NULL )
			x = x->left;
		return;
	}

	// look for the least ancestor where x was the left descendant
	while( x->parent != NULL){
		if( x->parent->left == x ){
			x = x->parent;
			return;
		}
		x = x->parent;
	}

	// x is already the right-most tree value
	x = NULL;
}

template< class Key, class Allocator >
void GisSplayTree< Key, Allocator >::decrement( GisNode*& x) {
	// if x is null find the rightmost node
	if( x == NULL ){
		x = root;
		while( x->right != NULL )
			x = x->right;
		return;
	}

	// if x has a left child, find its rightmost descendant
	if( x->left != NULL ){
		x = x->left;
		while( x->right != NULL )
			x = x->right;
		return;
	}

	// look for the least ancestor where x was the right descendant
	while( x->parent != NULL ){
		if( x->parent->right == x ){
			x = x->parent;
			return;
		}
		x = x->parent;
	}

	// x is already the left-most tree value
	// throw an exception?
	throw "decrementing out of bounds!";
	x = NULL;
}

template< class Key, class Allocator >
void GisSplayTree< Key, Allocator >::deleteTree() {
	if( root == NULL )
		return;
	// don't create a potentially huge call stack
	size_type zero = 0;
	GisNode* next_node = recursiveFind( zero, root );
	while( next_node != NULL ){
		GisNode* cur_node = next_node;
		increment( next_node );
		// unlink cur_node from the tree
		if( cur_node->right != NULL )
			cur_node->right->parent = cur_node->parent;
		if( cur_node->parent != NULL )
			cur_node->parent->left = cur_node->right;
		// sanity check
		if( cur_node->left != NULL )
			std::cerr << "insane in the brain\n";
		// delete cur_node
		delete cur_node->key;
		delete cur_node;
	}
}


#endif // __GisSplayTree_h__
