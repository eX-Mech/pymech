import wx
from wx import glcanvas
import OpenGL.GL as gl
from math import sqrt, atan2, asin, cos, sin


class MeshFrame(wx.Frame):
    """
    A frame to display meshes
    """

    def __init__(
        self,
        mesh,
        parent,
        id,
        title,
        pos=wx.DefaultPosition,
        size=wx.DefaultSize,
        style=wx.DEFAULT_FRAME_STYLE,
        name="frame",
    ):
        super(MeshFrame, self).__init__(parent, id, title, pos, size, style, name)
        self.GLinitialized = False
        attribList = (
            glcanvas.WX_GL_RGBA,  # RGBA
            glcanvas.WX_GL_DOUBLEBUFFER,  # Double Buffered
            glcanvas.WX_GL_DEPTH_SIZE,
            24,
        )  # 24 bit

        # Create the canvas
        self.canvas = glcanvas.GLCanvas(self, attribList=attribList)
        self.context = glcanvas.GLContext(self.canvas)

        # Set the event handlers.
        self.canvas.Bind(wx.EVT_ERASE_BACKGROUND, self.processEraseBackgroundEvent)
        self.canvas.Bind(wx.EVT_SIZE, self.processSizeEvent)
        self.canvas.Bind(wx.EVT_PAINT, self.processPaintEvent)

        # create a menu bar
        # self.makeMenuBar()

        # mesh display properties
        self.curve_points = 12

        # data to be drawn
        self.vertices = []
        self.edges = []
        self.buildMesh(mesh)

        # view parameters
        self.margins = 0.05
        # sets self.mesh_limits
        self.setLimits(mesh)
        # current limits
        self.limits = self.mesh_limits

        # and a status bar
        # self.CreateStatusBar()
        # self.SetStatusText("initialised")

    # Canvas Proxy Methods

    def GetGLExtents(self):
        """Get the extents of the OpenGL canvas."""
        return self.canvas.GetClientSize()

    def SwapBuffers(self):
        """Swap the OpenGL buffers."""
        self.canvas.SwapBuffers()

    def OnExit(self, event):
        """Close the frame, terminating the application."""
        self.Close(True)

    # wxPython Window Handlers

    def processEraseBackgroundEvent(self, event):
        """Process the erase background event."""
        pass  # Do nothing, to avoid flashing on MSWin

    def processSizeEvent(self, event):
        """Process the resize event."""
        if self.context:
            # Make sure the frame is shown before calling SetCurrent.
            self.Show()
            self.canvas.SetCurrent(self.context)

            size = self.GetGLExtents()
            self.OnReshape(size.width, size.height)
            self.canvas.Refresh(False)
        event.Skip()

    def processPaintEvent(self, event):
        """Process the drawing event."""
        self.canvas.SetCurrent(self.context)

        # This is a 'perfect' time to initialize OpenGL ... only if we need to
        if not self.GLinitialized:
            self.OnInitGL()
            self.GLinitialized = True

        self.OnDraw()
        event.Skip()

    # GLFrame OpenGL Event Handlers

    def OnInitGL(self):
        """Initialize OpenGL for use in the window."""
        gl.glClearColor(1, 1, 1, 1)

    def OnReshape(self, width, height):
        """Reshape the OpenGL viewport based on the dimensions of the window."""

        xmin = self.limits[0]
        xmax = self.limits[1]
        ymin = self.limits[2]
        ymax = self.limits[3]
        # check whether the view is limited by width or height, and scale accordingly
        lx = xmax - xmin
        ly = ymax - ymin
        if lx / width > ly / height:
            y0 = 0.5 * (ymin + ymax)
            dy = height / width * lx / 2
            ymin = y0 - dy
            ymax = y0 + dy
        else:
            x0 = 0.5 * (xmin + xmax)
            dx = width / height * ly / 2
            xmin = x0 - dx
            xmax = x0 + dx
        gl.glViewport(0, 0, width, height)
        gl.glMatrixMode(gl.GL_PROJECTION)
        gl.glLoadIdentity()
        gl.glOrtho(xmin, xmax, ymin, ymax, -1, 1)

        gl.glMatrixMode(gl.GL_MODELVIEW)
        gl.glLoadIdentity()

    def OnDraw(self, *args, **kwargs):
        "Draw the window."
        gl.glClear(gl.GL_COLOR_BUFFER_BIT)
        gl.glEnable(
            gl.GL_LINE_SMOOTH
        )  # this doesn't seem to be doing anything? It would be nice to have antialiasing
        gl.glLineWidth(1.0)
        gl.glBegin(gl.GL_LINES)
        gl.glColor(0, 0, 0)
        for edge in self.edges:
            for vertex in edge:
                gl.glVertex3fv(self.vertices[vertex])
        gl.glEnd()

        self.SwapBuffers()

    def buildMesh(self, mesh):
        """
        Builds the edges to be drawn based on the mesh representation.
        """

        # gives the indices of the vertices of an element in a position array
        def vertex_indices(iface):
            if iface == 0:
                return (0, 0)
            elif iface == 1:
                return (0, -1)
            elif iface == 2:
                return (-1, -1)
            else:
                return (-1, 0)

        current_point = 0
        first_point = 0
        for el in mesh.elem:
            first_point = current_point
            for iface in range(4):
                j0, i0 = vertex_indices(iface)
                if el.ccurv[iface] == "":
                    self.vertices.append(
                        (
                            el.pos[0, 0, j0, i0],
                            el.pos[1, 0, j0, i0],
                            0.0,
                        )
                    )
                    if iface < 3:
                        next_point = current_point + 1
                    else:
                        next_point = first_point
                    self.edges.append((current_point, next_point))
                    current_point += 1
                elif el.ccurv[iface] == "m":
                    # we should draw a parabola passing through the current vertex, the midpoint, and the next vertex.
                    x0, y0 = el.pos[0:2, 0, j0, i0]
                    xm, ym = el.curv[iface][0:2]
                    j1, i1 = vertex_indices((iface + 1) % 4)
                    x1, y1 = el.pos[0:2, 0, j1, i1]
                    # quadratic Lagrange interpolation between points
                    for ipt in range(self.curve_points):
                        # tp varies between 0 and 1
                        tp = ipt / self.curve_points
                        xp = (
                            x0 * 2 * (tp - 0.5) * (tp - 1)
                            - xm * 4 * tp * (tp - 1)
                            + x1 * 2 * tp * (tp - 0.5)
                        )
                        yp = (
                            y0 * 2 * (tp - 0.5) * (tp - 1)
                            - ym * 4 * tp * (tp - 1)
                            + y1 * 2 * tp * (tp - 0.5)
                        )
                        self.vertices.append((xp, yp, 0))
                        if iface == 3 and ipt == self.curve_points - 1:
                            next_point = first_point
                        else:
                            next_point = current_point + 1
                        self.edges.append((current_point, next_point))
                        current_point += 1
                elif el.ccurv[iface] == "C":
                    # draw a circle of given radius passing through the next vertex and the current one
                    # first, find the distance between the midpoint of the segment ((x0, y0), (x1, y1)) and the center (xc, yc) of the circle
                    radius = el.curv[iface][
                        0
                    ]  # this can be positive or negative depending on direction
                    x0, y0 = el.pos[0:2, 0, j0, i0]
                    j1, i1 = vertex_indices((iface + 1) % 4)
                    x1, y1 = el.pos[0:2, 0, j1, i1]
                    # length of the segment
                    ls2 = (x1 - x0) ** 2 + (y1 - y0) ** 2
                    try:
                        dist = radius * sqrt(1 - ls2 / (4 * radius**2))
                    except ValueError:
                        raise ValueError("the radius of the curved edge is too small")
                    # midpoint of the edge
                    xm = 0.5 * (x0 + x1)
                    ym = 0.5 * (y0 + y1)
                    # outward normal direction
                    ls = sqrt(ls2)
                    nx = (y1 - y0) / ls
                    ny = -(x1 - x0) / ls
                    # position of the centre
                    xc = xm - nx * dist
                    yc = ym - ny * dist
                    # now find the range of arguments spanned by the circle arc
                    # argument to the centre of the edge
                    theta0 = atan2(ym - yc, xm - xc)
                    dtheta = asin(ls / (2 * radius))
                    theta_min = theta0 - dtheta
                    # Now, add the points
                    for itheta in range(self.curve_points):
                        theta = theta_min + 2 * dtheta * itheta / self.curve_points
                        xp = xc + abs(radius) * cos(theta)
                        yp = yc + abs(radius) * sin(theta)
                        self.vertices.append((xp, yp, 0))
                        if iface == 3 and itheta == self.curve_points - 1:
                            next_point = first_point
                        else:
                            next_point = current_point + 1
                        self.edges.append((current_point, next_point))
                        current_point += 1

    def setLimits(self, mesh):
        """
        set view limits to the size of the mesh with some margin
        """
        xmin, xmax = mesh.lims.pos[0]
        ymin, ymax = mesh.lims.pos[1]
        lx = xmax - xmin
        ly = ymax - ymin
        self.mesh_limits = [
            xmin - self.margins * lx,
            xmax + self.margins * lx,
            ymin - self.margins * ly,
            ymax + self.margins * ly,
        ]


def plot2D(mesh):
    # make a new app & frame
    app = wx.App()
    frame = MeshFrame(mesh, None, -1, title="pymech")

    frame.Show()

    # Start the event loop.
    app.MainLoop()
