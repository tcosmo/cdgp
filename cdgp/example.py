""" This module is an example module with example function.
"""

def square(x: int) -> int:
  """ This function squares numbers. You can use it like that: 
  :code:`square(5)`.

    Args:
        `x` (int): number to square

    Returns:
        int: squared number

    :Example:
        >>> square(3)
        9
        >>> square(10)
        100
  """
  return x**2
