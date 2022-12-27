import wx
from wx import glcanvas
import OpenGL
from OpenGL.GL import *
from OpenGL.GLU import *

class MeshFrame(wx.Frame):
    """
    A frame to display meshes
    """

    def __init__(self,
                 mesh,
                 parent,
                 id,
                 title,
                 pos=wx.DefaultPosition,
                 size=wx.DefaultSize,
                 style=wx.DEFAULT_FRAME_STYLE,
                 name='frame'
    ):
        super(MeshFrame, self).__init__(parent, id, title, pos, size, style, name)
        self.GLinitialized = False
        attribList = (glcanvas.WX_GL_RGBA, # RGBA
                      glcanvas.WX_GL_DOUBLEBUFFER, # Double Buffered
                      glcanvas.WX_GL_DEPTH_SIZE, 24) # 24 bit

        # Create the canvas
        self.canvas = glcanvas.GLCanvas(self, attribList=attribList)
        self.context = glcanvas.GLContext(self.canvas)

        # Set the event handlers.
        self.canvas.Bind(wx.EVT_ERASE_BACKGROUND, self.processEraseBackgroundEvent)
        self.canvas.Bind(wx.EVT_SIZE, self.processSizeEvent)
        self.canvas.Bind(wx.EVT_PAINT, self.processPaintEvent)


        # create a menu bar
        #self.makeMenuBar()

        # data to be drawn
        self.vertices = []
        self.edges = []
        self.buildMesh(mesh)

        # view parameters
        self.margins = 0.05
        self.setLimits(mesh)  # sets xmin, xmax, ymin, ymax
        
        # and a status bar
        #self.CreateStatusBar()
        #self.SetStatusText("initialised")

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
        pass # Do nothing, to avoid flashing on MSWin

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
        glClearColor(1, 1, 1, 1)

    def OnReshape(self, width, height):
        """Reshape the OpenGL viewport based on the dimensions of the window."""
        glViewport(0, 0, width, height)

        glMatrixMode(GL_PROJECTION)
        glLoadIdentity()
        glOrtho(
            self.limits[0],
            self.limits[1],
            self.limits[2],
            self.limits[3],
            -1,
            1
        )

        glMatrixMode(GL_MODELVIEW)
        glLoadIdentity()

    def OnDraw(self, *args, **kwargs):
        "Draw the window."
        glClear(GL_COLOR_BUFFER_BIT)
        glEnable(GL_LINE_SMOOTH)  # this doesn't seem to be doing anything? It would be nice to have antialiasing
        glLineWidth(1.0)
        glBegin(GL_LINES)
        glColor(0, 0, 0)
        for edge in self.edges:
            for vertex in edge:
                glVertex3fv(self.vertices[vertex])
        glEnd()

        self.SwapBuffers()

    def buildMesh(self, mesh):
        k = 0
        for el in mesh.elem:
            self.vertices.append((
                el.pos[0, 0, 0, 0],
                el.pos[1, 0, 0, 0],
                0.,
            ))
            self.vertices.append((
                el.pos[0, 0, 0, -1],
                el.pos[1, 0, 0, -1],
                0.,
            ))
            self.vertices.append((
                el.pos[0, 0, -1, -1],
                el.pos[1, 0, -1, -1],
                0.,
            ))
            self.vertices.append((
                el.pos[0, 0, -1, 0],
                el.pos[1, 0, -1, 0],
                0.,
            ))
            self.edges.append((k, k + 1))
            self.edges.append((k + 1, k + 2))
            self.edges.append((k + 2, k + 3))
            self.edges.append((k + 3, k))
            k += 4

    def setLimits(self, mesh):
        """
        set view limits to the size of the mesh with some margin
        """
        xmin, xmax = mesh.lims.pos[0]
        ymin, ymax = mesh.lims.pos[1]
        lx = xmax - xmin
        ly = ymax - ymin
        self.limits = [
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
