import ipywidgets as ipw


class Spinner(ipw.VBox):
    """
    A spinner widget.
    """

    def __init__(self, *args, spinner_file=None, message="", **kwargs):
        self._spinner_file = spinner_file
        super().__init__(*args, **kwargs)
        with open(spinner_file) as f:
            self._spinner = ipw.HTML(f.read())
        self._message = ipw.HTML(message)
        self.children = [self._message, self._spinner]
        self.layout.display = "none"
        self.layout.width = "200px"

    @property
    def message(self):
        return self._message.value

    @property
    def spinner_file(self):
        return self._spinner_file

    def start(self):
        self.layout.display = "flex"

    def stop(self):
        self.layout.display = "none"
