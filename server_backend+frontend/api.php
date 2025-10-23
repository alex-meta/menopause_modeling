<?php
header('Content-Type: application/json; charset=utf-8');

$logDir = __DIR__ . '/log';
if (!is_dir($logDir)) { @mkdir($logDir, 0777, true); }

function log_write($name, $text) {
  global $logDir;
  $ts = date('Ymd_His');
  file_put_contents("$logDir/{$name}_$ts.log", $text);
}

try {
  if (!isset($_POST['payload'])) {
    echo json_encode(['error' => 'payload missing']); exit;
  }

  $payload = $_POST['payload'];
  $test = json_decode($payload, true);
  if ($test === null) { echo json_encode(['error' => 'invalid JSON']); exit; }

  $R = '/usr/bin/Rscript';                         // проверь путь
  $script = __DIR__ . '/swan_predictor2.R';        // проверь путь

  if (!file_exists($script)) {
    echo json_encode(['error' => 'R script not found']); exit;
  }

  $arg = escapeshellarg($payload);
  $cmd = escapeshellcmd($R) . ' --vanilla ' . escapeshellarg($script) . ' ' . $arg . ' 2>&1';

  log_write('request', "CMD: $cmd\nPAYLOAD:\n$payload\n");

  $output = shell_exec($cmd);
  if ($output === null) {
    echo json_encode(['error' => 'failed to execute Rscript']); exit;
  }

  log_write('r_output', $output);

  $json = json_decode($output, true);
  if ($json === null) {
    echo json_encode(['error' => 'R output is not JSON', 'raw' => $output]); exit;
  }

  echo json_encode($json);
} catch (Throwable $e) {
  log_write('php_error', $e->getMessage() . "\n" . $e->getTraceAsString());
  echo json_encode(['error' => $e->getMessage()]);
}
