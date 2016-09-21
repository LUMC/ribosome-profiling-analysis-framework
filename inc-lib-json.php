<?php
/*******************************************************************************
 *
 * This script provides the JSON related functions for PHP versions that do not
 * contain these, i.e. PHP < 5.2.0.
 *
 *
 * Code taken with permission from:
 * https://gajendrakrjain.wordpress.com/2014/10/07/alternative-for-json_decode-and-json_encode/
 *
 * This work is licensed under the Creative Commons
 * Attribution-NonCommercial-ShareAlike 4.0 International License. To view a
 * copy of this license, visit http://creativecommons.org/licenses/by-nc-sa/4.0/
 * or send a letter to:
 * Creative Commons, PO Box 1866, Mountain View, CA 94042, USA.
 *
 *************/

if (!function_exists('json_encode')) {
    function json_encode($a = false)
    {
        // Some basic debugging to ensure we have something returned.
        if (is_null($a)) return 'null';
        if ($a === false) return 'false';
        if ($a === true) return 'true';
        if (is_scalar($a)) {
            if (is_float($a)) {
                // Always use '.' for floats.
                return floatval(str_replace(',', '.', strval($a)));
            }
            if (is_string($a)) {
                static $jsonReplaces = array(array('\\', '/', "\n", "\t", "\r", "\b", "\f", '"'), array('\\\\', '\\/', '\\n', '\\t', '\\r', '\\b', '\\f', '\"'));
                return '"' . str_replace($jsonReplaces[0], $jsonReplaces[1], $a) . '"';
            } else
                return $a;
        }
        $isList = true;
        for ($i = 0, reset($a); true; $i++) {
            if (key($a) !== $i) {
                $isList = false;
                break;
            }
        }
        $result = array();
        if ($isList) {
            foreach ($a as $v) $result[] = json_encode($v);
            return '[' . join(',', $result) . ']';
        } else {
            foreach ($a as $k => $v) $result[] = json_encode($k) . ':' . json_encode($v);
            return '{' . join(',', $result) . '}';
        }
    }
}

if (!function_exists('json_decode')) {
    function json_decode($json)
    {
        $comment = false;
        $out = '$x=';
        for ($i = 0; $i < strlen($json); $i++) {
            if (!$comment) {
                if (($json[$i] == '{') || ($json[$i] == '[')) $out .= ' array('; else if (($json[$i] == '}') || ($json[$i] == ']')) $out .= ')'; else if ($json[$i] == ':') $out .= '=>';
                else
                    $out .= $json[$i];
            } else
                $out .= $json[$i];
            if ($json[$i] == '"' && $json[($i - 1)] != "\\")
                $comment = !$comment;
        }
        eval($out . ';');
        return $x;
    }
}
?>
