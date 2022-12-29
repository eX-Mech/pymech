import wx
from wx import glcanvas
import OpenGL.GL as gl
import numpy as np
from math import sqrt, atan2, asin, cos, sin
import time


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
        self.canvas.Bind(wx.EVT_MOUSEWHEEL, self.processMouseWheelEvent)

        # create a menu bar
        # self.makeMenuBar()

        # mesh display properties
        self.curve_points = 12

        # data to be drawn
        self.vertex_data = np.array([])
        self.colour_data = np.array([])
        self.num_vertices = 0
        self.buildMesh(mesh)
        self.vertex_buffer = 0
        self.colour_buffer = 0

        # view parameters
        # relative margins to display around the mesh in the default view
        self.margins = 0.05
        # when zooming in by one step, the field of view is multiplied by this value
        self.zoom_factor = 0.95
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
            self.updateLimits(size.width, size.height)
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

    def processMouseWheelEvent(self, event):
        """
        Processes mouse wheel events.
        Zooms when the vertical mouse wheel is used.
        """

        # vertical or horizontal wheel axis
        axis = event.GetWheelAxis()
        # logical context for the mouse position
        dc = wx.ClientDC(self)
        # position of the pointer in pixels in the window frame
        xcp, ycp = event.GetLogicalPosition(dc).Get()
        increment = event.GetWheelRotation() / event.GetWheelDelta()
        if axis == wx.MOUSE_WHEEL_VERTICAL:
            self.zoom(xcp, ycp, increment)

    def zoom(self, xcp, ycp, increment):
        """
        Zooms around the given centre proportionally to the increment.
        Zooms in if `increment` is positive, out if it is negative.
        The centre (xcp, ycp) is given in pixels, not in mesh units.
        """

        factor = self.zoom_factor**increment
        xmin, xmax, ymin, ymax = self.limits
        size = self.GetGLExtents()
        # get the centre location in mesh units
        # xcp is 0 on the left, ycp is zero on the top
        xc = xmin + (xmax - xmin) * xcp / size.width
        yc = ymax - (ymax - ymin) * ycp / size.height
        # update the limits
        xmin = xc + factor * (xmin - xc)
        xmax = xc + factor * (xmax - xc)
        ymin = yc + factor * (ymin - yc)
        ymax = yc + factor * (ymax - yc)
        self.limits = [xmin, xmax, ymin, ymax]
        size = self.GetGLExtents()
        self.updateLimits(size.width, size.height)
        self.OnDraw()

    # GLFrame OpenGL Event Handlers

    def OnInitGL(self):
        """Initialize OpenGL for use in the window."""
        self.createBuffers()
        gl.glClearColor(1, 1, 1, 1)

    def updateLimits(self, width, height):
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

    def createBuffers(self):
        # new vertex buffer
        self.vertex_buffer = gl.glGenBuffers(1)
        gl.glBindBuffer(gl.GL_ARRAY_BUFFER, self.vertex_buffer)
        # send the vertex data to the buffer
        gl.glBufferData(gl.GL_ARRAY_BUFFER, self.vertex_data, gl.GL_STATIC_DRAW)
        # unbind the buffer
        gl.glBindBuffer(gl.GL_ARRAY_BUFFER, 0)

    def OnDraw(self, *args, **kwargs):
        "Draw the window."

        t1 = time.perf_counter()
        # initialise
        gl.glClear(gl.GL_COLOR_BUFFER_BIT)
        gl.glClear(gl.GL_DEPTH_BUFFER_BIT)
        gl.glClearColor(1, 1, 1, 1)
        gl.glEnable(gl.GL_LINE_SMOOTH)
        gl.glLineWidth(1.0)
        gl.glEnableClientState(gl.GL_VERTEX_ARRAY)
        gl.glColor(0, 0, 0)
        # load buffers
        gl.glBindBuffer(gl.GL_ARRAY_BUFFER, self.vertex_buffer)
        gl.glVertexPointer(3, gl.GL_DOUBLE, 0, None)
        gl.glBindBuffer(gl.GL_ARRAY_BUFFER, 0)
        # draw the mesh
        gl.glDrawArrays(gl.GL_LINES, 0, self.num_vertices)
        # finalise
        gl.glDisableClientState(gl.GL_VERTEX_ARRAY)
        t2 = time.perf_counter()

        self.SwapBuffers()
        print(f"time: draw {t2 - t1:.6e}")

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
        vertices = []
        edges = []
        for el in mesh.elem:
            first_point = current_point
            for iface in range(4):
                j0, i0 = vertex_indices(iface)
                if el.ccurv[iface] == "":
                    vertices.append(
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
                    edges.append((current_point, next_point))
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
                        vertices.append((xp, yp, 0))
                        if iface == 3 and ipt == self.curve_points - 1:
                            next_point = first_point
                        else:
                            next_point = current_point + 1
                        edges.append((current_point, next_point))
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
                        vertices.append((xp, yp, 0))
                        if iface == 3 and itheta == self.curve_points - 1:
                            next_point = first_point
                        else:
                            next_point = current_point + 1
                        edges.append((current_point, next_point))
                        current_point += 1

        # put everything into a buffer that OpenGL can read
        self.num_vertices = 2 * len(edges)
        self.colour_data = np.array([0 for _ in range(4 * self.num_vertices)])
        vertex_data = []
        for edge in edges:
            for vertex in edge:
                x, y, z = vertices[vertex]
                vertex_data.append(x)
                vertex_data.append(y)
                vertex_data.append(z)
        self.vertex_data = np.array(vertex_data)

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
