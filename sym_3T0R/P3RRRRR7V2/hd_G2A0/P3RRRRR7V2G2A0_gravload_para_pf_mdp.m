% Calculate minimal parameter regressor of Gravitation load for parallel robot
% P3RRRRR7V2G2A0
% Use Code from Maple symbolic Code Generation
%
% Input:
% xP [3x1]
%   Generalized platform coordinates
% qJ [3x3]
%   Generalized joint coordinates (joint angles)
%   rows: links of the robot
%   columns: number of leg
% legFrame [3x3]
%   base frame orientation for each leg
%   row: number of leg
%   column: Euler angles for the orientation.
%   Euler angle convention from robot definition ("leg_frame")
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,d3,d4]';
% koppelP [3x3]
%   coordinates of the platform coupling joints
%   (joints that link the end of legs with platform)
%   in platform coordinates
%   rows: number of leg
%   columns: x-, y-, z-coordinates
% MDP [18x1]
%   Minimal dynamic parameter vector for parallel robot(fixed base model)
%   see P3RRRRR7V2G2A0_convert_par2_MPV_fixb.m

% Output:
% taugX [3x1]
%   minimal parameter regressor of Gravitation load for parallel robot
%   in task space

% Quelle: HybrDyn-Toolbox
% Datum: 2020-08-07 10:08
% Revision: 8f4ff0ee124033641e65b154ac60823cef59ef1f (2020-07-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taugX = P3RRRRR7V2G2A0_gravload_para_pf_mdp(xP, qJ, g, legFrame, ...
  koppelP, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,3),zeros(3,1),zeros(3,3),zeros(3,3),zeros(7,1),zeros(18,1)}
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3RRRRR7V2G2A0_gravload_para_pf_mdp: xP has to be [3x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 3]), ...
  'P3RRRRR7V2G2A0_gravload_para_pf_mdp: qJ has to be [3x3] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'P3RRRRR7V2G2A0_gravload_para_pf_mdp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'P3RRRRR7V2G2A0_gravload_para_pf_mdp: pkin has to be [7x1] (double)');
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3RRRRR7V2G2A0_gravload_para_pf_mdp: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3RRRRR7V2G2A0_gravload_para_pf_mdp: Koppelpunkt has to be [3x3] (double)');
assert(isreal(MDP) && all(size(MDP) == [18 1]), ...
  'P3RRRRR7V2G2A0_gravload_para_pf_mdp: MDP has to be [18x1] (double)'); 

%% Symbolic Calculation
% From invdyn_para_plfcoord_taugreg_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2020-08-07 10:08:54
% EndTime: 2020-08-07 10:08:59
% DurationCPUTime: 5.01s
% Computational Cost: add. (2196->360), mult. (3711->541), div. (210->14), fcn. (2796->66), ass. (0->266)
t8402 = legFrame(1,2);
t8366 = sin(t8402);
t8369 = cos(t8402);
t8336 = t8369 * g(1) - t8366 * g(2);
t8411 = sin(qJ(1,1));
t8420 = cos(qJ(1,1));
t8312 = -g(3) * t8411 + t8336 * t8420;
t8397 = qJ(2,1) + qJ(3,1);
t8357 = cos(t8397);
t8419 = cos(qJ(2,1));
t8330 = 0.1e1 / (t8419 * pkin(2) + pkin(3) * t8357 + pkin(1));
t8528 = t8330 * t8411;
t8474 = t8312 * t8528;
t8575 = 0.2e1 * pkin(3);
t8385 = (pkin(5) + pkin(6) + pkin(7));
t8574 = 2 * t8385;
t8413 = cos(qJ(2,3));
t8390 = t8413 ^ 2;
t8573 = 0.2e1 * t8390;
t8416 = cos(qJ(2,2));
t8392 = t8416 ^ 2;
t8572 = 0.2e1 * t8392;
t8394 = t8419 ^ 2;
t8571 = 0.2e1 * t8394;
t8412 = cos(qJ(3,3));
t8389 = t8412 ^ 2;
t8378 = pkin(3) * t8389;
t8415 = cos(qJ(3,2));
t8391 = t8415 ^ 2;
t8379 = pkin(3) * t8391;
t8418 = cos(qJ(3,1));
t8393 = t8418 ^ 2;
t8380 = pkin(3) * t8393;
t8403 = sin(qJ(3,3));
t8569 = t8403 * pkin(1);
t8568 = t8403 * pkin(3);
t8406 = sin(qJ(3,2));
t8567 = t8406 * pkin(1);
t8566 = t8406 * pkin(3);
t8409 = sin(qJ(3,1));
t8565 = t8409 * pkin(1);
t8564 = t8409 * pkin(3);
t8563 = t8412 * pkin(2);
t8370 = t8412 * pkin(3);
t8562 = t8415 * pkin(2);
t8371 = t8415 * pkin(3);
t8561 = t8418 * pkin(2);
t8372 = t8418 * pkin(3);
t8560 = -qJ(3,1) + qJ(1,1);
t8559 = qJ(3,1) + qJ(1,1);
t8558 = -qJ(3,2) + qJ(1,2);
t8557 = qJ(3,2) + qJ(1,2);
t8556 = -qJ(3,3) + qJ(1,3);
t8555 = qJ(3,3) + qJ(1,3);
t8554 = qJ(1,1) + 0.2e1 * qJ(2,1);
t8553 = qJ(1,1) - 0.2e1 * qJ(2,1);
t8552 = qJ(1,2) + 0.2e1 * qJ(2,2);
t8551 = qJ(1,2) - 0.2e1 * qJ(2,2);
t8550 = 0.2e1 * qJ(2,3) + qJ(1,3);
t8549 = -0.2e1 * qJ(2,3) + qJ(1,3);
t8548 = 0.2e1 * pkin(1);
t8395 = qJ(2,3) + qJ(3,3);
t8352 = sin(t8395);
t8358 = qJ(2,3) + t8555;
t8359 = -qJ(2,3) + t8556;
t8404 = sin(qJ(2,3));
t8414 = cos(qJ(1,3));
t8424 = 0.2e1 * qJ(3,3);
t8435 = pkin(2) ^ 2;
t8547 = ((cos(t8359) + cos(t8358)) * t8548 + (sin(t8359) + sin(t8358)) * t8574 + (cos(0.2e1 * qJ(3,3) - t8549) + cos(t8424 + t8550) + 0.2e1 * t8414) * pkin(3) + (cos(qJ(3,3) - t8549) + cos(qJ(3,3) + t8550) + cos(t8556) + cos(t8555)) * pkin(2)) / (-t8435 * sin(qJ(2,3) - qJ(3,3)) + pkin(2) * (0.2e1 * t8569 + pkin(2) * t8352 + (sin(t8424 + qJ(2,3)) - t8404) * pkin(3)));
t8396 = qJ(2,2) + qJ(3,2);
t8353 = sin(t8396);
t8360 = qJ(2,2) + t8557;
t8361 = -qJ(2,2) + t8558;
t8407 = sin(qJ(2,2));
t8417 = cos(qJ(1,2));
t8427 = 0.2e1 * qJ(3,2);
t8546 = ((cos(t8361) + cos(t8360)) * t8548 + (sin(t8361) + sin(t8360)) * t8574 + (cos(0.2e1 * qJ(3,2) - t8551) + cos(t8427 + t8552) + 0.2e1 * t8417) * pkin(3) + (cos(qJ(3,2) - t8551) + cos(qJ(3,2) + t8552) + cos(t8558) + cos(t8557)) * pkin(2)) / (-t8435 * sin(qJ(2,2) - qJ(3,2)) + pkin(2) * (0.2e1 * t8567 + pkin(2) * t8353 + (sin(t8427 + qJ(2,2)) - t8407) * pkin(3)));
t8354 = sin(t8397);
t8362 = qJ(2,1) + t8559;
t8363 = -qJ(2,1) + t8560;
t8410 = sin(qJ(2,1));
t8430 = 0.2e1 * qJ(3,1);
t8545 = ((cos(t8363) + cos(t8362)) * t8548 + (sin(t8363) + sin(t8362)) * t8574 + (cos(0.2e1 * qJ(3,1) - t8553) + cos(t8430 + t8554) + 0.2e1 * t8420) * pkin(3) + (cos(qJ(3,1) - t8553) + cos(qJ(3,1) + t8554) + cos(t8560) + cos(t8559)) * pkin(2)) / (-t8435 * sin(qJ(2,1) - qJ(3,1)) + pkin(2) * (0.2e1 * t8565 + pkin(2) * t8354 + (sin(t8430 + qJ(2,1)) - t8410) * pkin(3)));
t8405 = sin(qJ(1,3));
t8454 = pkin(1) * t8405 - t8414 * t8385;
t8504 = t8403 * t8404;
t8298 = t8454 * t8504 + (t8389 - 0.1e1) * t8405 * pkin(3);
t8400 = legFrame(3,2);
t8364 = sin(t8400);
t8544 = t8298 * t8364;
t8367 = cos(t8400);
t8543 = t8298 * t8367;
t8408 = sin(qJ(1,2));
t8453 = pkin(1) * t8408 - t8417 * t8385;
t8503 = t8406 * t8407;
t8299 = t8453 * t8503 + (t8391 - 0.1e1) * t8408 * pkin(3);
t8401 = legFrame(2,2);
t8365 = sin(t8401);
t8542 = t8299 * t8365;
t8368 = cos(t8401);
t8541 = t8299 * t8368;
t8452 = pkin(1) * t8411 - t8420 * t8385;
t8502 = t8409 * t8410;
t8300 = t8452 * t8502 + (t8393 - 0.1e1) * t8411 * pkin(3);
t8540 = t8300 * t8366;
t8539 = t8300 * t8369;
t8334 = t8367 * g(1) - t8364 * g(2);
t8316 = -g(3) * t8405 + t8334 * t8414;
t8355 = cos(t8395);
t8328 = 0.1e1 / (t8413 * pkin(2) + pkin(3) * t8355 + pkin(1));
t8538 = t8316 * t8328;
t8335 = t8368 * g(1) - t8365 * g(2);
t8317 = -g(3) * t8408 + t8335 * t8417;
t8356 = cos(t8396);
t8329 = 0.1e1 / (t8416 * pkin(2) + pkin(3) * t8356 + pkin(1));
t8537 = t8317 * t8329;
t8536 = t8312 * t8330;
t8491 = pkin(3) * t8504;
t8349 = t8370 + pkin(2);
t8520 = t8349 * t8413;
t8535 = 0.1e1 / (pkin(1) - t8491 + t8520) / t8403;
t8490 = pkin(3) * t8503;
t8350 = t8371 + pkin(2);
t8517 = t8350 * t8416;
t8534 = 0.1e1 / (pkin(1) - t8490 + t8517) / t8406;
t8489 = pkin(3) * t8502;
t8351 = t8372 + pkin(2);
t8514 = t8351 * t8419;
t8533 = 0.1e1 / (pkin(1) - t8489 + t8514) / t8409;
t8532 = t8328 * t8405;
t8531 = t8328 * t8414;
t8530 = t8329 * t8408;
t8529 = t8329 * t8417;
t8527 = t8330 * t8420;
t8423 = pkin(2) / 0.2e1;
t8525 = (t8370 + t8423) * t8403;
t8524 = (t8371 + t8423) * t8406;
t8523 = (t8372 + t8423) * t8409;
t8522 = t8349 * t8364;
t8521 = t8349 * t8367;
t8519 = t8350 * t8365;
t8518 = t8350 * t8368;
t8516 = t8351 * t8366;
t8515 = t8351 * t8369;
t8513 = t8364 * t8405;
t8512 = t8365 * t8408;
t8511 = t8366 * t8411;
t8510 = t8367 * t8405;
t8509 = t8368 * t8408;
t8508 = t8369 * t8411;
t8433 = pkin(3) ^ 2;
t8507 = t8389 * t8433;
t8506 = t8391 * t8433;
t8505 = t8393 * t8433;
t8343 = pkin(1) * t8404 - t8568;
t8501 = t8412 * t8343;
t8344 = pkin(1) * t8407 - t8566;
t8500 = t8415 * t8344;
t8345 = pkin(1) * t8410 - t8564;
t8499 = t8418 * t8345;
t8498 = -t8433 / 0.2e1 + t8435 / 0.2e1;
t8497 = pkin(2) * t8370;
t8496 = pkin(2) * t8371;
t8495 = pkin(2) * t8372;
t8494 = t8349 * t8568;
t8493 = t8350 * t8566;
t8492 = t8351 * t8564;
t8331 = t8364 * g(1) + t8367 * g(2);
t8457 = g(3) * t8414 + t8334 * t8405;
t8283 = -t8331 * t8355 + t8352 * t8457;
t8488 = t8283 * t8535;
t8284 = t8331 * t8352 + t8355 * t8457;
t8487 = t8284 * t8535;
t8332 = t8365 * g(1) + t8368 * g(2);
t8456 = g(3) * t8417 + t8335 * t8408;
t8285 = -t8332 * t8356 + t8353 * t8456;
t8486 = t8285 * t8534;
t8286 = t8332 * t8353 + t8356 * t8456;
t8485 = t8286 * t8534;
t8333 = t8366 * g(1) + t8369 * g(2);
t8455 = g(3) * t8420 + t8336 * t8411;
t8287 = -t8333 * t8357 + t8354 * t8455;
t8484 = t8287 * t8533;
t8288 = t8333 * t8354 + t8357 * t8455;
t8483 = t8288 * t8533;
t8289 = -t8331 * t8413 + t8404 * t8457;
t8482 = t8289 * t8535;
t8290 = t8331 * t8404 + t8413 * t8457;
t8481 = t8290 * t8535;
t8291 = -t8332 * t8416 + t8407 * t8456;
t8480 = t8291 * t8534;
t8292 = t8332 * t8407 + t8416 * t8456;
t8479 = t8292 * t8534;
t8293 = -t8333 * t8419 + t8410 * t8455;
t8478 = t8293 * t8533;
t8294 = t8333 * t8410 + t8419 * t8455;
t8477 = t8294 * t8533;
t8476 = t8413 * t8538;
t8475 = t8416 * t8537;
t8473 = t8312 * t8527;
t8472 = t8355 * t8538;
t8471 = t8364 * t8531;
t8470 = t8367 * t8531;
t8469 = t8356 * t8537;
t8468 = t8365 * t8529;
t8467 = t8368 * t8529;
t8466 = t8354 * t8536;
t8465 = t8357 * t8536;
t8464 = t8366 * t8527;
t8463 = t8369 * t8527;
t8462 = t8316 * t8532;
t8461 = t8405 * t8504;
t8460 = t8317 * t8530;
t8459 = t8408 * t8503;
t8458 = t8411 * t8502;
t8451 = t8414 * t8476;
t8450 = t8417 * t8475;
t8449 = t8410 * t8473;
t8448 = t8419 * t8473;
t8447 = t8316 * t8471;
t8446 = t8317 * t8468;
t8445 = t8414 * t8472;
t8444 = t8316 * t8470;
t8443 = t8417 * t8469;
t8442 = t8317 * t8467;
t8441 = t8420 * t8466;
t8440 = t8420 * t8465;
t8313 = t8461 * t8575 - t8454;
t8439 = pkin(2) * t8461 + t8313 * t8412;
t8314 = t8459 * t8575 - t8453;
t8438 = pkin(2) * t8459 + t8314 * t8415;
t8315 = t8458 * t8575 - t8452;
t8437 = pkin(2) * t8458 + t8315 * t8418;
t8436 = 0.1e1 / pkin(2);
t8434 = 0.1e1 / pkin(3);
t8422 = -pkin(3) / 0.2e1;
t8384 = -t8433 + t8435;
t8339 = t8380 + t8561 / 0.2e1 + t8422;
t8338 = t8379 + t8562 / 0.2e1 + t8422;
t8337 = t8378 + t8563 / 0.2e1 + t8422;
t8327 = t8495 + t8498 + t8505;
t8326 = t8496 + t8498 + t8506;
t8325 = t8497 + t8498 + t8507;
t8306 = t8565 + (-pkin(3) + t8561 + 0.2e1 * t8380) * t8410;
t8305 = t8567 + (-pkin(3) + t8562 + 0.2e1 * t8379) * t8407;
t8304 = t8569 + (-pkin(3) + t8563 + 0.2e1 * t8378) * t8404;
t8303 = pkin(1) * t8564 + (t8384 + 0.2e1 * t8495 + 0.2e1 * t8505) * t8410;
t8302 = pkin(1) * t8566 + (t8384 + 0.2e1 * t8496 + 0.2e1 * t8506) * t8407;
t8301 = pkin(1) * t8568 + (t8384 + 0.2e1 * t8497 + 0.2e1 * t8507) * t8404;
t8282 = -0.2e1 * t8327 * t8420 * t8394 - ((pkin(1) - 0.2e1 * t8489) * t8420 + t8411 * t8385) * t8514 + pkin(3) * ((pkin(1) * t8502 - pkin(3) + t8380) * t8420 + t8385 * t8458);
t8281 = -0.2e1 * t8326 * t8417 * t8392 - ((pkin(1) - 0.2e1 * t8490) * t8417 + t8408 * t8385) * t8517 + pkin(3) * ((pkin(1) * t8503 - pkin(3) + t8379) * t8417 + t8385 * t8459);
t8280 = -0.2e1 * t8325 * t8414 * t8390 - ((pkin(1) - 0.2e1 * t8491) * t8414 + t8405 * t8385) * t8520 + pkin(3) * ((pkin(1) * t8504 - pkin(3) + t8378) * t8414 + t8385 * t8461);
t8276 = (t8339 * t8508 + t8366 * t8523) * t8571 + (t8366 * t8306 - t8437 * t8369) * t8419 - t8539 + t8366 * t8499;
t8275 = (-t8339 * t8511 + t8369 * t8523) * t8571 + (t8369 * t8306 + t8437 * t8366) * t8419 + t8540 + t8369 * t8499;
t8274 = (t8338 * t8509 + t8365 * t8524) * t8572 + (t8365 * t8305 - t8438 * t8368) * t8416 - t8541 + t8365 * t8500;
t8273 = (-t8338 * t8512 + t8368 * t8524) * t8572 + (t8368 * t8305 + t8438 * t8365) * t8416 + t8542 + t8368 * t8500;
t8272 = (t8337 * t8510 + t8364 * t8525) * t8573 + (t8364 * t8304 - t8439 * t8367) * t8413 - t8543 + t8364 * t8501;
t8271 = (-t8337 * t8513 + t8367 * t8525) * t8573 + (t8367 * t8304 + t8439 * t8364) * t8413 + t8544 + t8367 * t8501;
t8270 = (-t8327 * t8508 - t8366 * t8492) * t8571 + (-t8366 * t8303 + t8315 * t8515) * t8419 + pkin(3) * t8539 - t8345 * t8516;
t8269 = (t8327 * t8511 - t8369 * t8492) * t8571 + (-t8303 * t8369 - t8315 * t8516) * t8419 - pkin(3) * t8540 - t8345 * t8515;
t8268 = (-t8326 * t8509 - t8365 * t8493) * t8572 + (-t8365 * t8302 + t8314 * t8518) * t8416 + pkin(3) * t8541 - t8344 * t8519;
t8267 = (t8326 * t8512 - t8368 * t8493) * t8572 + (-t8302 * t8368 - t8314 * t8519) * t8416 - pkin(3) * t8542 - t8344 * t8518;
t8266 = (-t8325 * t8510 - t8364 * t8494) * t8573 + (-t8364 * t8301 + t8313 * t8521) * t8413 + pkin(3) * t8543 - t8343 * t8522;
t8265 = (t8325 * t8513 - t8367 * t8494) * t8573 + (-t8301 * t8367 - t8313 * t8522) * t8413 - pkin(3) * t8544 - t8343 * t8521;
t1 = [(-t8312 * t8463 - t8442 - t8444) * MDP(2) + (t8455 * t8463 + t8456 * t8467 + t8457 * t8470) * MDP(3) + (-t8367 * t8451 - t8368 * t8450 - t8369 * t8448) * MDP(9) + (t8369 * t8449 + t8404 * t8444 + t8407 * t8442) * MDP(10) + (-t8367 * t8445 - t8368 * t8443 - t8369 * t8440) * MDP(16) + (t8352 * t8444 + t8353 * t8442 + t8369 * t8441) * MDP(17) - g(1) * MDP(18) + ((t8272 * t8482 + t8274 * t8480 + t8276 * t8478) * MDP(9) + (t8272 * t8481 + t8274 * t8479 + t8276 * t8477) * MDP(10) + (t8272 * t8488 + t8274 * t8486 + t8276 * t8484) * MDP(16) + (t8272 * t8487 + t8274 * t8485 + t8276 * t8483) * MDP(17) + ((t8266 * t8488 + t8268 * t8486 + t8270 * t8484) * MDP(16) + (t8266 * t8487 + t8268 * t8485 + t8270 * t8483) * MDP(17)) * t8434) * t8436; (t8312 * t8464 + t8446 + t8447) * MDP(2) + (-t8455 * t8464 - t8456 * t8468 - t8457 * t8471) * MDP(3) + (t8364 * t8451 + t8365 * t8450 + t8366 * t8448) * MDP(9) + (-t8366 * t8449 - t8404 * t8447 - t8407 * t8446) * MDP(10) + (t8364 * t8445 + t8365 * t8443 + t8366 * t8440) * MDP(16) + (-t8352 * t8447 - t8353 * t8446 - t8366 * t8441) * MDP(17) - g(2) * MDP(18) + ((t8271 * t8482 + t8273 * t8480 + t8275 * t8478) * MDP(9) + (t8271 * t8481 + t8273 * t8479 + t8275 * t8477) * MDP(10) + (t8271 * t8488 + t8273 * t8486 + t8275 * t8484) * MDP(16) + (t8271 * t8487 + t8273 * t8485 + t8275 * t8483) * MDP(17) + ((t8265 * t8488 + t8267 * t8486 + t8269 * t8484) * MDP(16) + (t8265 * t8487 + t8267 * t8485 + t8269 * t8483) * MDP(17)) * t8434) * t8436; (t8460 + t8462 + t8474) * MDP(2) + (-t8455 * t8528 - t8456 * t8530 - t8457 * t8532) * MDP(3) + (t8405 * t8476 + t8408 * t8475 + t8419 * t8474) * MDP(9) + (-t8404 * t8462 - t8407 * t8460 - t8410 * t8474) * MDP(10) + (t8405 * t8472 + t8408 * t8469 + t8411 * t8465) * MDP(16) + (-t8352 * t8462 - t8353 * t8460 - t8411 * t8466) * MDP(17) - g(3) * MDP(18) + ((t8280 * t8488 + t8281 * t8486 + t8282 * t8484) * MDP(16) + (t8280 * t8487 + t8281 * t8485 + t8282 * t8483) * MDP(17)) * t8436 * t8434 + (t8289 * t8547 + t8291 * t8546 + t8293 * t8545) * MDP(9) / 0.2e1 + (t8290 * t8547 + t8292 * t8546 + t8294 * t8545) * MDP(10) / 0.2e1 + (t8283 * t8547 + t8285 * t8546 + t8287 * t8545) * MDP(16) / 0.2e1 + (t8284 * t8547 + t8286 * t8546 + t8288 * t8545) * MDP(17) / 0.2e1;];
taugX  = t1;
