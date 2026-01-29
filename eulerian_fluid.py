import math
import random
import tkinter
import numpy

#some presets:
#x = 10, y = 6, box = 100, arrow = 30, font = 16
#x = 30, y = 15, box = 40, arrow = 40, font = 8
#x = 50, y = 30, box = 20, arrow = 30, font = 5
xCells = 30
yCells = 15
boxSize = 40
arrowSize = 40
fontSize = 8

fracWall = 0.2
inPressure = 100
outPressure = 0
deltaTime = 0.02
density = 1000
length = 0.01
viscosity = 0.1

grid = numpy.full( (xCells, yCells), True)
xVelocityGrid = numpy.full( (xCells + 1, yCells), True)
yVelocityGrid = numpy.full( (xCells, yCells + 1), True)
xVelocities = numpy.zeros( (xCells + 1, yCells))
yVelocities = numpy.zeros( (xCells, yCells + 1))
pressures = numpy.zeros( (xCells, yCells))

xArrows = numpy.empty( (xCells + 1, yCells), dtype = object)
yArrows = numpy.empty( (xCells, yCells + 1), dtype = object)
pressureDisplays = numpy.empty( (xCells, yCells), dtype = object)
pressureColors = numpy.empty( (xCells, yCells), dtype = object)

root = tkinter.Tk()
root.geometry(f"{xCells * boxSize + 100}x{yCells * boxSize + 100}")
root.title("Lattice Boltzmann Flow")
canvas = tkinter.Canvas(root, bg = "gray6")
canvas.pack(fill = "both", expand = True)

#updates grid, xVelocityGrid, yVelocityGrid to make walls and solid cells
def createGrid():
    #making walls on top and bottom
    for i in range(xCells):
        yVelocityGrid[i][0] = False
        yVelocityGrid[i][yCells] = False
    #filling in solid cells
    for i in range(xCells):
        for j in range(yCells):
            if random.random() < fracWall:
                grid[i][j] = False
                xVelocityGrid[i][j] = False
                xVelocityGrid[i + 1][j] = False
                yVelocityGrid[i][j] = False
                yVelocityGrid[i][j + 1] = False
    for i in range(xCells):
        for j in range(yCells):
            pressures[i][j] = (random.random()-.5)*10
                
#returns the pressure of a cell, or the inlet or outlet pressure
def getPressure(xCoord, yCoord):
    if xCoord == -1:
        return inPressure
    if xCoord == xCells:
        return outPressure
    return pressures[xCoord][yCoord]

#returns edge count, sum of neighbouring pressures, and divergence a cell
def getPressureCalculationVariables(xCoord, yCoord):
    edgeCount = 0
    pressureSum = 0
    divergence = 0
    if xVelocityGrid[xCoord][yCoord]:
        edgeCount += 1
        pressureSum += getPressure(xCoord - 1, yCoord)
        divergence -= xVelocities[xCoord][yCoord]
    if xVelocityGrid[xCoord + 1][yCoord]:
        edgeCount += 1
        pressureSum += getPressure(xCoord + 1, yCoord)
        divergence += xVelocities[xCoord + 1][yCoord]
    if yVelocityGrid[xCoord][yCoord]:
        edgeCount += 1
        pressureSum += getPressure(xCoord, yCoord - 1)
        divergence -= yVelocities[xCoord][yCoord]
    if yVelocityGrid[xCoord][yCoord + 1]:
        edgeCount += 1
        pressureSum += getPressure(xCoord, yCoord + 1)
        divergence += yVelocities[xCoord][yCoord + 1]
    return [edgeCount, pressureSum, divergence]

#updates the velocity according the surrounding pressures and the flow in/out of the cell
def updatePressures():
    temp = numpy.zeros( (xCells, yCells) )
    for i in range(xCells):
        for j in range(yCells):
            [edgeCount, pressureSum, divergence] = getPressureCalculationVariables(i, j)
            if edgeCount == 0:
                temp[i][j] = 0
            else:
                temp[i][j] = (pressureSum - density * length * divergence / deltaTime) / edgeCount
    global pressures
    pressures = temp

#finds laplacian of x velocities
def getXVelocityLaplacian(xCoord, yCoord):
    velocity = xVelocities[xCoord][yCoord]
    nearbySum = 0
    if yCoord > 0 and xVelocityGrid[xCoord][yCoord - 1]:
        nearbySum += xVelocities[xCoord][yCoord - 1]
    if yCoord < yCells - 1 and xVelocityGrid[xCoord][yCoord + 1]:
        nearbySum += xVelocities[xCoord][yCoord + 1]
    if xCoord > 0 and xVelocityGrid[xCoord - 1][yCoord]:
        nearbySum += xVelocities[xCoord - 1][yCoord]
    if xCoord < xCells and xVelocityGrid[xCoord + 1][yCoord]:
        nearbySum += xVelocities[xCoord + 1][yCoord]
    if xCoord == 0 or xCoord == xCells:
        nearbySum += velocity
    return (nearbySum - 4 * velocity) / (length ** 2)

#finds laplacian of y velocities
def getYVelocityLaplacian(xCoord, yCoord):
    velocity = yVelocities[xCoord][yCoord]
    nearbySum = 0
    if xCoord > 0 and yVelocityGrid[xCoord - 1][yCoord]:
        nearbySum += yVelocities[xCoord - 1][yCoord]
    if xCoord < xCells - 1 and yVelocityGrid[xCoord + 1][yCoord]:
        nearbySum += yVelocities[xCoord + 1][yCoord]
    if yCoord > 0 and yVelocityGrid[xCoord][yCoord - 1]:
        nearbySum += yVelocities[xCoord][yCoord - 1]
    if yCoord < yCells and yVelocityGrid[xCoord][yCoord + 1]:
        nearbySum += yVelocities[xCoord][yCoord + 1]
    return (nearbySum - 4 * velocity) / (length ** 2)

#updates velocities depending on pressure drop and viscosity
def updateVelocities():
    xDeltaV = numpy.zeros( (xCells + 1, yCells) )
    yDeltaV = numpy.zeros( (xCells, yCells + 1) )
    for i in range(xCells + 1):
        for j in range(yCells):
            if xVelocityGrid[i][j]:
                leftPressure = getPressure(i - 1, j)
                rightPressure = getPressure(i, j)
                laplacian = getXVelocityLaplacian(i, j)
                xDeltaV[i][j] = ((leftPressure - rightPressure) / length  + viscosity * laplacian) * deltaTime / density
                
    for i in range(xCells):
        for j in range(yCells + 1):
            if yVelocityGrid[i][j]:
                topPressure = getPressure(i, j - 1)
                bottomPressure = getPressure(i, j)
                laplacian = getYVelocityLaplacian(i, j)
                yDeltaV[i][j] = ((topPressure - bottomPressure) / length  + viscosity * laplacian) * deltaTime / density
    global xVelocities
    global yVelocities
    xVelocities += xDeltaV
    yVelocities += yDeltaV

#draws the grid, pressure displays, and velocity arrows
def drawGrid(): 
    for i in range(xCells):
        for j in range(yCells):
            if grid[i][j]:
                pressureColors[i][j] = canvas.create_rectangle(i * boxSize + 50, j * boxSize + 50, (i + 1) * boxSize + 50, (j + 1) * boxSize + 50, outline = "", fill = "white")
                pressureDisplays[i][j] = canvas.create_text((i + 0.5) * boxSize + 50, (j + 0.5) * boxSize + 50, text = f"{pressures[i][j]: .2f}", font=("Helvetica", fontSize), fill = "gray6")
            else:
                canvas.create_rectangle(i * boxSize + 50, j * boxSize + 50, (i + 1) * boxSize + 50, (j + 1) * boxSize + 50, outline = "", fill = "gray16")
    
    canvas.create_rectangle(49, 49, xCells * boxSize + 51, yCells * boxSize + 51, outline = "white", width = 2)
    for i in range(1, xCells):
        canvas.create_line(i * boxSize + 50, 50, i * boxSize + 50, yCells * boxSize + 50, fill = "gray16")
    for i in range(1, yCells):
        canvas.create_line(50, i * boxSize + 50, xCells * boxSize + 50, i * boxSize + 50, fill = "gray16")

    #currently drawing arrows into solid cells, could change
    for i in range(xCells + 1):
        for j in range(yCells):
            xArrows[i][j] = canvas.create_line(i * boxSize + 50, (j + 0.5) * boxSize + 50, i * boxSize + 50 + xVelocities[i][j] * arrowSize, (j + 0.5) * boxSize + 50, arrow = tkinter.LAST, arrowshape = (5, 3, 2),fill = "SpringGreen3", width = 2)
    for i in range(xCells):
        for j in range(yCells + 1):
            yArrows[i][j] = canvas.create_line((i + 0.5) * boxSize + 50, j * boxSize + 50, (i + 0.5) * boxSize + 50, j * boxSize + 50 + yVelocities[i][j] * arrowSize, arrow = tkinter.LAST, arrowshape = (5, 3, 2),fill = "SpringGreen3", width = 2)

#updates pressure display, color, and velocity arrows
def updateDisplay():
    for i in range(xCells):
        for j in range(yCells):
            if grid[i][j]:
                pressure = pressures[i][j]
                percent = (pressure - outPressure) / (inPressure - outPressure)
                canvas.itemconfig(pressureDisplays[i][j], text = f"{pressure: .2f}")
                if percent <= 0:
                    t = (percent + 1.5) / 1.5
                    r = 255
                    g = int(255 * t)
                    b = int(255 * t)
                else:
                    t = percent / 1.5
                    r = int(255 * (1 - t))
                    g = int(255 * (1 - t))
                    b = 255
                canvas.itemconfig(pressureColors[i][j], fill = f"#{r:02x}{g:02x}{b:02x}")
    
    #currently drawing arrows into solid cells, could change
    for i in range(xCells + 1):
        for j in range(yCells):
            canvas.coords(xArrows[i][j], i * boxSize + 50, (j + 0.5) * boxSize + 50, i * boxSize + 50 + xVelocities[i][j] * arrowSize, (j + 0.5) * boxSize + 50)
    for i in range(xCells):
        for j in range(yCells + 1):
            canvas.coords(yArrows[i][j], (i + 0.5) * boxSize + 50, j * boxSize + 50, (i + 0.5) * boxSize + 50, j * boxSize + 50 + yVelocities[i][j] * arrowSize)

def runSimulation():
    for i in range(10):
        updatePressures()
    updateVelocities()
    updateDisplay()
    root.after(20, runSimulation)

createGrid()
drawGrid()
runSimulation()
root.mainloop()