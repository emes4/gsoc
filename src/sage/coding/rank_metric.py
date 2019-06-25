from sage.coding.linear_code import AbstractLinearCode
from sage.rings.integer import Integer
from sage.coding.relative_finite_field_extension import *
#import sage.coding.gabidulin

def to_matrix_representation(base_field, sub_field, v):
    """
    """
    if not v.is_vector():
        raise TypeError("Input must be a vector")
    n = v.length()
    FE = RelativeFiniteFieldExtension(Fqm, Fq)
    m = Fqm.degree()//Fq.degree()
    g = matrix(Fq, m, n, lambda i,j: FE.relative_field_representation(v[j])[i])
    return g

def from_matrix_representation(base_field, sub_field, m):
    """
    """
    if not m.is_matrix():
        raise TypeError("Input must be a matrix")
    FE = RelativeFiniteFieldExtension(Fqm, Fq)
    v = []
    for i in range(m.ncols()):
        v.append(FE.absolute_field_representation(m.column(i)))
    return vector(v)

def rank_weight(base_field, sub_field, c):
    """
    """
    if c.is_vector(c):
        c = _to_matrix_representation(Fqm, Fq, c)
    return c.rank()

def rank_distance(base_field, sub_field, a, b):
    """
    """
    if a.is_vector():
        a = _to_matrix_representation(Fqm, Fq, a)
    if b.is_vector():
        b = _to_matrix_representation(Fqm, Fq, b)
    return (a - b).rank()


class AbstractLinearRankMetricCode(AbstractLinearCode):
    """
    Abstract class for linear rank metric codes.

    This class contains methods that can be used on families of linear rank
    metric codes. Every linear rank metric code class should inherit from this
    abstract class.

    For details on how to implement a linear rank metric code, follow the
    instructions in documentation for
    :class:`sage.coding.linear_code.AbstractLinearCode`.
    (Or should I just copy and paste the text from there?)
    """

    _registered_encoders = {}
    _registered_decoders = {}

    def __init__(self, base_field, sub_field, length, dimension, \
            default_encoder_name, default_decoder_name, field_extension=None):
        """
        Initializes mandatory parameters that every linear rank metric code has.

        This method only exists for inheritance purposes as it initializes
        parameters that need to be known by every linear rank metric code.
        The class :class:`sage.coding.rank_metric.AbstractLinearRankMetricCode`
        should never be directly instantiated.

        INPUT:

        - ``base_field`` -- the base field of ``self``

        - ``sub_field`` -- the sub field of ``self``

        - ``length`` -- the length of ``self`` (a Python int or a Sage Integer,
          must be > 0)

        - ``dimension`` -- the dimension of ``self``

        - ``field_extension`` -- representation of the elements of the relative
          extension of `base_field` over `sub_field` (default: ``None``)

        - ``default_encoder_name`` -- the name of the default encoder of ``self``

        - ``default_decoder_name`` -- the name of the default decoder of ``self``
        """
        if not isinstance(dimension, (int, Integer)):
            raise ValueError("dimension must be a Python int or a Sage Integer")
        if not sub_field.is_field():
            raise ValueError("'sub_field' must be a field (and {} is not one)".format(sub_field))
        if not field_extension: #if field_extension is provided, then what? how to check?
            field_extension = RelativeFiniteFieldExtension(base_field, sub_field)
        self._base_field = base_field
        self._sub_field = sub_field
        self._length = length
        self._dimension = dimension
        self._field_extension = field_extension

        #TODO: super()

    def base_field(self):
        """
        Returns the base field of ``self``.
        """
        return self._base_field

    def sub_field(self):
        """
        Returns the sub field of ``self``.
        """
        return self._sub_field

    def field_extension(self):
        """
        Returns the field extension of ``self``.
        """
        return self._field_extension

    def distance(self, left, right):
        """
        Returns the rank of the matrix of ``left`` - ``right``.
        """
        return rank_distance(self._base_field, self._sub_field, left, right)

    def minimum_distance(self):
        r"""
        Return an error requiring to override ``minimum_distance`` in ``self``.

        There is currently no general algorithm calculating the minimum distance
        of linear rank metric codes. One has to implement the specific method
        when writing a new code class which inherits from
        :class:`AbstractLinearRankMetricCode`.
        The generic call to ``minimum_distance`` has to fail.
        """

        #TODO: check that the formatting self.parent() works
        raise RuntimeError("Please override minimum_distance in the implementation of {}".format(self.parent()))

    def weight(self, word):
        """
        Returns the weight of the code word - its rank.
        """
        return rank_weight(self._base_field, self._sub_field, word)

    def to_matrix(self, word):
        """
        Returns the matrix representation of a word.
        """
        return to_matrix_representation(self._base_field, self._sub_field, word)

    def from_matrix(self, word):
        """
        Returns the vector representation of a word.
        """
        return from_matrix_representation(self._base_field, self._sub_field, word)

    @property
    def automorphism_group_gens(self):
        raise AttributeError("%s has no attribute 'automorphism_group_gens'" % self.__class__)

    @property
    def canonical_representative(self):
        raise AttributeError("%s has no attribute 'canonical_representative'" % self.__class__)

    @property
    def permutation_automorphism_group(self):
        raise AttributeError("%s has no attribute 'permutation_automorphism_group'" % self.__class__)

    @property
    def assmuss_mattson_designs(self):
        raise AttributeError("%s has no attribute 'assmus_mattson_designs'" % self.__class__)

    @property
    def binomial_moment(self):
        raise AttributeError("%s has no attribute 'binomial_moment'" % self.__class__)
