def multiply_array(arr):
    result = 1
    for num in arr:
        result *= num
    return result

# Example usage
numbers = [7, 7, 73, 127, 337, 92737, 649657]
product = multiply_array(numbers)
print(f"The product of {numbers} is: {product}")

# Alternative shorter version using reduce