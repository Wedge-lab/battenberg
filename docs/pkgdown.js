/* http://gregfranko.com/blog/jquery-best-practices/ */
(($) => {
  $(() => {
    $(".navbar-fixed-top").headroom();

    const updateBodyPadding = () => {
      $("body").css("padding-top", $(".navbar").height() + 10);
    };

    updateBodyPadding();
    $(window).resize(updateBodyPadding);

    $('[data-toggle="tooltip"]').tooltip();

    const cur_path = paths(location.pathname);
    const links = $("#navbar ul li a");
    let max_length = -1;
    let pos = -1;

    links.each((i, link) => {
      if (link.getAttribute("href") === "#") return;
      if (link.host !== location.host) return;

      const nav_path = paths(link.pathname);
      const length = prefix_length(nav_path, cur_path);

      if (length > max_length) {
        max_length = length;
        pos = i;
      }
    });

    if (pos >= 0) {
      const menu_anchor = $(links[pos]);
      menu_anchor.parent().addClass("active");
      menu_anchor.closest("li.dropdown").addClass("active");
    }
  });

  const paths = (pathname) => {
    const pieces = pathname.split("/");
    pieces.shift(); // always starts with /

    const end = pieces[pieces.length - 1];
    if (end === "index.html" || end === "") pieces.pop();
    return pieces;
  };

  const prefix_length = (needle, haystack) => {
    if (needle.length > haystack.length) return -1;
    if (haystack.length === 0) return needle.length === 0 ? 0 : -1;

    for (let i = 0; i < haystack.length; i++) {
      if (needle[i] !== haystack[i]) return i;
    }
    return haystack.length;
  };

  /* Clipboard --------------------------*/

  const changeTooltipMessage = (element, msg) => {
    const tooltipOriginalTitle = element.getAttribute("data-original-title");
    element.setAttribute("data-original-title", msg);
    $(element).tooltip("show");
    element.setAttribute("data-original-title", tooltipOriginalTitle);
  };

  if (window.ClipboardJS && ClipboardJS.isSupported()) {
    $(document).ready(() => {
      const copyButton =
        "<button type='button' class='btn btn-primary btn-copy-ex' title='Copy to clipboard' aria-label='Copy to clipboard' data-toggle='tooltip' data-placement='left auto' data-trigger='hover' data-clipboard-copy><i class='fa fa-copy'></i></button>";

      $("div.sourceCode").addClass("hasCopyButton").prepend(copyButton);

      $(".btn-copy-ex").tooltip({ container: "body" });

      const clipboardBtnCopies = new ClipboardJS("[data-clipboard-copy]", {
        text: (trigger) =>
          trigger.parentNode.textContent.replace(/\n#>[^\n]*/g, ""),
      });

      clipboardBtnCopies.on("success", (e) => {
        changeTooltipMessage(e.trigger, "Copied!");
        e.clearSelection();
      });

      clipboardBtnCopies.on("error", (e) => {
        changeTooltipMessage(e.trigger, "Press Ctrl+C or Command+C to copy");
      });
    });
  }
})(window.jQuery || window.$);
