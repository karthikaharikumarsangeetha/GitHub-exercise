import matplotlib.pyplot as plt

# Data
mortality_types = ['Overall Mortality', 'Liver-specific Mortality']
rates = [25.56, 11.77]

# Create bar chart
plt.figure(figsize=(8, 5))
bars = plt.bar(mortality_types, rates, color=['skyblue', 'salmon'])
plt.ylabel('Deaths per 1,000 patient-years')
plt.title('Mortality in Patients with Confirmed NASH')
plt.ylim(0, 30)

# Add data labels
for bar in bars:
    yval = bar.get_height()
    plt.text(bar.get_x() + bar.get_width()/2, yval + 0.5,
             f'{yval:.2f}', ha='center', fontweight='bold')

plt.show()
