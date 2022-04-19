/** Amy Chan
 * You may wish to copy the contents of this into your prism.min.js and minify it.
 * 27-05-2013: creation.
 *  New tokens you need to define in CSS: .token.(slot|namespace|function|variable)
 *  From `?Quotes`:
 *   Identifiers consist of a sequence of letters, digits, the period
 *    (‘.’) and the underscore.  They must not start with a digit nor
 *    underscore, nor with a period followed by a digit.  Reserved words
 *    are not valid identifiers.
 * 13-12-2013: modified.
 *  Fixed bug where namespaces wouldn't highlight.
 *  Included partial workaround for comments-within-strings bug (commented out
 *   by default).
 */
Prism.languages.r = {
	'comment': /#.*$/gm,    
	// if you want to have a partial workaround for the comments-within-strings bug, use the below instead of the above.
	// However a comment with an odd number of " in it will not highlight.
	//'comment': /#(?=(?:[^"\\\r\n]*(\\.|"(?:[^"\\\r\n]*\\.)*[^"\\\r\n]*"))*[^"\r\n]*$).*$/gm,
    'string': /("|')(?:\\.|(?!\\|\1)[\s\S])*\1/g,
    'keyword': /\b(?:if|else|repeat|while|function|for|in|next|break)\b/g,
    // NULL etc are not really booleans but I just group them tohere to be marked up
    'boolean': /\b(?:TRUE|FALSE|T|F|NA(?:_(?:integer|real|complex|character)_)?|NULL)\b/g,
    'function': /(?:(?:[a-zA-Z]|\.(?![0-9]))[.\w]*|`[^`\s]*`)[ \t]*(?=\()/g,
    'number': /\b[-+]?(0x[\dA-Fa-f]+|\d*\.?\d+([Ee]-?\d+)?i?|Inf|NaN)\b/g,
    'operator': /(?:<|&lt;)-|[-+\/!\^]|={1,2}|(?:&[lg]t;|[><])=?|(&amp;|&){1,2}|\|\|?|\*\*?|\%(\/|in)?\%/g,
    'property': {
        pattern: /([$@])[\w._]+/g,
        lookbehind: true
    },
	'namespace': /[\w._]+(?=::)/g,
    'punctuation': /[@?${}[\];(),.:]/g, // Prism.Languages.clike.punctuation
    'variable': /(?:(?:[a-zA-Z]|\.(?![0-9]))[.\w]*|`[^`\s]*`)/g
};