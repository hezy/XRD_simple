using Plots
gr()  # Initialize the GR backend

# Generate some data
x = 1:10
y = rand(10)

# Create the plot
p = plot(x, y, title="Random Data", xlabel="X-axis", ylabel="Y-axis")

# Display the plot
display(p)
