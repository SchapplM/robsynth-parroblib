% Calculate minimal parameter regressor of Gravitation load for parallel robot
% P4RRRRR10V2G1A0
% Use Code from Maple symbolic Code Generation
%
% Input:
% xP [4x1]
%   Generalized platform coordinates
% qJ [3x4]
%   Generalized joint coordinates (joint angles)
%   rows: links of the robot
%   columns: number of leg
% legFrame [4x3]
%   base frame orientation for each leg
%   row: number of leg
%   column: Euler angles for the orientation.
%   Euler angle convention from robot definition ("leg_frame")
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,alpha2,d1,d2,d3,d4]';
% koppelP [4x3]
%   coordinates of the platform coupling joints
%   (joints that link the end of legs with platform)
%   in platform coordinates
%   rows: number of leg
%   columns: x-, y-, z-coordinates

% Output:
% tau_reg [4x21]
%   minimal parameter regressor of Gravitation load for parallel robot
%   in task space

% Quelle: HybrDyn-Toolbox
% Datum: 2020-08-07 17:14
% Revision: 8f4ff0ee124033641e65b154ac60823cef59ef1f (2020-07-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = P4RRRRR10V2G1A0_gravload_para_pf_regmin(xP, qJ, g, legFrame, ...
  koppelP, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,4),zeros(3,1),zeros(4,3),zeros(4,3),zeros(8,1)}
assert(isreal(xP) && all(size(xP) == [4 1]), ...
  'P4RRRRR10V2G1A0_gravload_para_pf_regmin: xP has to be [4x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 4]), ...
  'P4RRRRR10V2G1A0_gravload_para_pf_regmin: qJ has to be [3x4] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'P4RRRRR10V2G1A0_gravload_para_pf_regmin: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'P4RRRRR10V2G1A0_gravload_para_pf_regmin: pkin has to be [8x1] (double)');
assert(isreal(legFrame) && all(size(legFrame) == [4 3]), ...
  'P4RRRRR10V2G1A0_gravload_para_pf_regmin: legFrame has to be [4x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [4 3]), ...
  'P4RRRRR10V2G1A0_gravload_para_pf_regmin: Koppelpunkt has to be [4x3] (double)');

%% Symbolic Calculation
% From invdyn_para_plfcoord_taugreg_matlab.m
% OptimizationMode: 2
% StartTime: 2020-08-07 15:39:44
% EndTime: 2020-08-07 15:40:24
% DurationCPUTime: 43.15s
% Computational Cost: add. (11393->1012), mult. (21608->1729), div. (240->29), fcn. (18920->36), ass. (0->692)
t8372 = cos(qJ(2,1));
t8374 = pkin(8) + pkin(7);
t8292 = t8374 * t8372;
t8363 = sin(qJ(2,1));
t8216 = pkin(2) * t8363 - t8292;
t8345 = cos(pkin(4));
t8702 = pkin(1) * t8345;
t8169 = t8216 * t8702;
t8344 = sin(pkin(4));
t8371 = cos(qJ(3,1));
t8289 = t8374 * t8363;
t8328 = t8372 * pkin(2);
t8219 = t8328 + t8289;
t8286 = t8328 + pkin(1);
t8362 = sin(qJ(3,1));
t8322 = t8362 * pkin(3);
t8407 = -pkin(6) * t8219 + t8286 * t8322;
t8201 = pkin(1) + t8219;
t8684 = t8201 * pkin(2);
t8492 = t8362 * t8684;
t8570 = t8345 * t8363;
t8581 = t8344 * t8372;
t8341 = t8371 ^ 2;
t8692 = pkin(3) * t8341;
t8423 = -(pkin(1) * t8570 + pkin(6) * t8581) * t8692 + t8344 * t8492;
t8735 = 0.1e1 / ((-t8344 * t8407 + t8169) * t8371 - t8423);
t8369 = cos(qJ(2,2));
t8291 = t8374 * t8369;
t8360 = sin(qJ(2,2));
t8215 = pkin(2) * t8360 - t8291;
t8168 = t8215 * t8702;
t8368 = cos(qJ(3,2));
t8288 = t8374 * t8360;
t8326 = t8369 * pkin(2);
t8218 = t8326 + t8288;
t8283 = t8326 + pkin(1);
t8359 = sin(qJ(3,2));
t8321 = t8359 * pkin(3);
t8408 = -pkin(6) * t8218 + t8283 * t8321;
t8200 = pkin(1) + t8218;
t8685 = t8200 * pkin(2);
t8494 = t8359 * t8685;
t8572 = t8345 * t8360;
t8584 = t8344 * t8369;
t8339 = t8368 ^ 2;
t8693 = pkin(3) * t8339;
t8424 = -(pkin(1) * t8572 + pkin(6) * t8584) * t8693 + t8344 * t8494;
t8734 = 0.1e1 / ((-t8344 * t8408 + t8168) * t8368 - t8424);
t8366 = cos(qJ(2,3));
t8290 = t8374 * t8366;
t8357 = sin(qJ(2,3));
t8214 = pkin(2) * t8357 - t8290;
t8167 = t8214 * t8702;
t8365 = cos(qJ(3,3));
t8287 = t8374 * t8357;
t8324 = t8366 * pkin(2);
t8217 = t8324 + t8287;
t8280 = t8324 + pkin(1);
t8356 = sin(qJ(3,3));
t8320 = t8356 * pkin(3);
t8409 = -pkin(6) * t8217 + t8280 * t8320;
t8199 = pkin(1) + t8217;
t8686 = t8199 * pkin(2);
t8496 = t8356 * t8686;
t8574 = t8345 * t8357;
t8587 = t8344 * t8366;
t8337 = t8365 ^ 2;
t8694 = pkin(3) * t8337;
t8425 = -(pkin(1) * t8574 + pkin(6) * t8587) * t8694 + t8344 * t8496;
t8733 = 0.1e1 / ((-t8344 * t8409 + t8167) * t8365 - t8425);
t8354 = cos(qJ(2,4));
t8271 = t8374 * t8354;
t8351 = sin(qJ(2,4));
t8212 = pkin(2) * t8351 - t8271;
t8165 = t8212 * t8702;
t8353 = cos(qJ(3,4));
t8270 = t8374 * t8351;
t8319 = t8354 * pkin(2);
t8213 = t8319 + t8270;
t8262 = t8319 + pkin(1);
t8350 = sin(qJ(3,4));
t8317 = t8350 * pkin(3);
t8410 = -pkin(6) * t8213 + t8262 * t8317;
t8198 = pkin(1) + t8213;
t8687 = t8198 * pkin(2);
t8498 = t8350 * t8687;
t8577 = t8345 * t8351;
t8596 = t8344 * t8354;
t8333 = t8353 ^ 2;
t8695 = pkin(3) * t8333;
t8426 = -(pkin(1) * t8577 + pkin(6) * t8596) * t8695 + t8344 * t8498;
t8732 = 0.1e1 / ((-t8344 * t8410 + t8165) * t8353 - t8426);
t8487 = pkin(6) * t8362 + pkin(3);
t8719 = -t8487 + t8692;
t8488 = pkin(6) * t8359 + pkin(3);
t8718 = -t8488 + t8693;
t8489 = pkin(6) * t8356 + pkin(3);
t8717 = -t8489 + t8694;
t8490 = pkin(6) * t8350 + pkin(3);
t8716 = -t8490 + t8695;
t8318 = t8353 * pkin(3);
t8260 = t8318 + pkin(2);
t8323 = t8365 * pkin(3);
t8278 = t8323 + pkin(2);
t8325 = t8368 * pkin(3);
t8281 = t8325 + pkin(2);
t8327 = t8371 * pkin(3);
t8284 = t8327 + pkin(2);
t8385 = 0.1e1 / pkin(3);
t8579 = t8344 * t8385;
t8258 = t8317 + pkin(6);
t8334 = t8354 ^ 2;
t8566 = t8345 * t8374;
t8507 = 0.2e1 * t8566;
t8431 = (t8334 - 0.1e1 / 0.2e1) * t8507;
t8335 = pkin(2) + t8374;
t8336 = pkin(2) - t8374;
t8604 = t8335 * t8336;
t8454 = t8345 * t8604;
t8462 = t8258 * t8596;
t8711 = t8213 * t8258 * t8344 - t8351 * t8354 * t8454 + pkin(2) * t8431 + (t8431 + t8462) * t8318;
t8272 = t8320 + pkin(6);
t8338 = t8366 ^ 2;
t8294 = t8338 - 0.1e1 / 0.2e1;
t8707 = -0.2e1 * pkin(2);
t8436 = t8566 * t8707;
t8460 = t8272 * t8587;
t8710 = t8272 * t8217 * t8344 - t8357 * t8366 * t8454 - t8294 * t8436 + (t8294 * t8507 + t8460) * t8323;
t8274 = t8321 + pkin(6);
t8340 = t8369 ^ 2;
t8295 = t8340 - 0.1e1 / 0.2e1;
t8459 = t8274 * t8584;
t8709 = t8274 * t8218 * t8344 - t8360 * t8369 * t8454 - t8295 * t8436 + (t8295 * t8507 + t8459) * t8325;
t8276 = t8322 + pkin(6);
t8342 = t8372 ^ 2;
t8296 = t8342 - 0.1e1 / 0.2e1;
t8458 = t8276 * t8581;
t8708 = t8219 * t8276 * t8344 - t8363 * t8372 * t8454 - t8296 * t8436 + (t8296 * t8507 + t8458) * t8327;
t8706 = pkin(1) * t8260;
t8705 = pkin(1) * t8278;
t8704 = pkin(1) * t8281;
t8703 = pkin(1) * t8284;
t8701 = pkin(2) * t8344;
t8700 = pkin(2) * t8345;
t8699 = pkin(2) * t8350;
t8698 = pkin(2) * t8356;
t8697 = pkin(2) * t8359;
t8696 = pkin(2) * t8362;
t8691 = pkin(6) * t8353;
t8690 = pkin(6) * t8365;
t8689 = pkin(6) * t8368;
t8688 = pkin(6) * t8371;
t8332 = t8345 ^ 2;
t8683 = (t8332 - 0.1e1) * pkin(6);
t8682 = pkin(3) * t8707;
t8352 = sin(qJ(1,4));
t8355 = cos(qJ(1,4));
t8376 = koppelP(4,2);
t8380 = koppelP(4,1);
t8204 = t8352 * t8380 - t8355 * t8376;
t8205 = t8352 * t8376 + t8355 * t8380;
t8375 = xP(4);
t8330 = sin(t8375);
t8331 = cos(t8375);
t8088 = t8204 * t8330 + t8205 * t8331;
t8346 = legFrame(4,3);
t8307 = sin(t8346);
t8311 = cos(t8346);
t8414 = t8204 * t8331 - t8205 * t8330;
t7999 = t8088 * t8311 - t8307 * t8414;
t8576 = t8345 * t8353;
t8113 = pkin(1) * (-t8260 * t8351 + t8271) * t8576;
t8173 = pkin(2) * t8334 + t8351 * t8271 - pkin(2);
t8301 = t8334 - 0.2e1;
t8550 = t8353 * t8354;
t8440 = pkin(3) * t8550 + t8213;
t8553 = t8351 * t8353;
t8599 = t8344 * t8350;
t8402 = -(t8262 * t8351 + t8374 + (pkin(3) * t8553 - t8271) * t8354) * t8599 + (pkin(1) + t8440) * t8576;
t8497 = t8351 * t8695;
t8513 = -t8212 * t8353 - 0.2e1 * t8497;
t8600 = t8344 * t8345;
t8681 = (t8402 * (t8088 * t8307 + t8311 * t8414) + (-t8513 * t8332 + (-t8490 * t8332 - t8716) * t8351 + ((t8301 * t8318 + t8173) * t8350 - t8691) * t8600) * t7999) / (t8113 + t8344 * ((t8262 * t8318 + t8687) * t8350 - t8440 * t8691));
t8361 = sin(qJ(1,2));
t8370 = cos(qJ(1,2));
t8378 = koppelP(2,2);
t8382 = koppelP(2,1);
t8208 = t8361 * t8382 - t8370 * t8378;
t8209 = t8361 * t8378 + t8370 * t8382;
t8096 = t8208 * t8330 + t8209 * t8331;
t8348 = legFrame(2,3);
t8309 = sin(t8348);
t8313 = cos(t8348);
t8412 = t8208 * t8331 - t8209 * t8330;
t8006 = t8096 * t8313 - t8309 * t8412;
t8568 = t8345 * t8368;
t8115 = pkin(1) * (-t8281 * t8360 + t8291) * t8568;
t8178 = pkin(2) * t8340 + t8360 * t8291 - pkin(2);
t8303 = t8340 - 0.2e1;
t8525 = t8368 * t8369;
t8438 = pkin(3) * t8525 + t8218;
t8537 = t8360 * t8368;
t8592 = t8344 * t8359;
t8400 = -(t8283 * t8360 + t8374 + (pkin(3) * t8537 - t8291) * t8369) * t8592 + (pkin(1) + t8438) * t8568;
t8493 = t8360 * t8693;
t8511 = -t8215 * t8368 - 0.2e1 * t8493;
t8680 = (t8400 * (t8309 * t8096 + t8313 * t8412) + (-t8511 * t8332 + (-t8332 * t8488 - t8718) * t8360 + ((t8303 * t8325 + t8178) * t8359 - t8689) * t8600) * t8006) / (t8115 + t8344 * ((t8283 * t8325 + t8685) * t8359 - t8438 * t8689));
t8364 = sin(qJ(1,1));
t8373 = cos(qJ(1,1));
t8379 = koppelP(1,2);
t8383 = koppelP(1,1);
t8210 = t8364 * t8383 - t8373 * t8379;
t8211 = t8364 * t8379 + t8373 * t8383;
t8097 = t8210 * t8330 + t8211 * t8331;
t8349 = legFrame(1,3);
t8310 = sin(t8349);
t8314 = cos(t8349);
t8411 = t8210 * t8331 - t8211 * t8330;
t8008 = t8097 * t8314 - t8310 * t8411;
t8567 = t8345 * t8371;
t8116 = pkin(1) * (-t8284 * t8363 + t8292) * t8567;
t8179 = pkin(2) * t8342 + t8363 * t8292 - pkin(2);
t8304 = t8342 - 0.2e1;
t8523 = t8371 * t8372;
t8437 = pkin(3) * t8523 + t8219;
t8530 = t8363 * t8371;
t8590 = t8344 * t8362;
t8399 = -(t8286 * t8363 + t8374 + (pkin(3) * t8530 - t8292) * t8372) * t8590 + (pkin(1) + t8437) * t8567;
t8491 = t8363 * t8692;
t8510 = -t8216 * t8371 - 0.2e1 * t8491;
t8679 = (t8399 * (t8310 * t8097 + t8314 * t8411) + (-t8510 * t8332 + (-t8332 * t8487 - t8719) * t8363 + ((t8304 * t8327 + t8179) * t8362 - t8688) * t8600) * t8008) / (t8116 + ((t8286 * t8327 + t8684) * t8362 - t8437 * t8688) * t8344);
t8358 = sin(qJ(1,3));
t8367 = cos(qJ(1,3));
t8377 = koppelP(3,2);
t8381 = koppelP(3,1);
t8206 = t8358 * t8381 - t8367 * t8377;
t8207 = t8358 * t8377 + t8367 * t8381;
t8095 = t8206 * t8330 + t8207 * t8331;
t8347 = legFrame(3,3);
t8308 = sin(t8347);
t8312 = cos(t8347);
t8413 = t8206 * t8331 - t8207 * t8330;
t8004 = t8095 * t8312 - t8308 * t8413;
t8569 = t8345 * t8365;
t8114 = pkin(1) * (-t8278 * t8357 + t8290) * t8569;
t8177 = pkin(2) * t8338 + t8357 * t8290 - pkin(2);
t8302 = t8338 - 0.2e1;
t8527 = t8365 * t8366;
t8439 = pkin(3) * t8527 + t8217;
t8544 = t8357 * t8365;
t8594 = t8344 * t8356;
t8401 = -(t8280 * t8357 + t8374 + (pkin(3) * t8544 - t8290) * t8366) * t8594 + (pkin(1) + t8439) * t8569;
t8495 = t8357 * t8694;
t8512 = -t8214 * t8365 - 0.2e1 * t8495;
t8678 = (t8401 * (t8095 * t8308 + t8312 * t8413) + (-t8512 * t8332 + (-t8489 * t8332 - t8717) * t8357 + ((t8302 * t8323 + t8177) * t8356 - t8690) * t8600) * t8004) / (t8114 + t8344 * ((t8280 * t8323 + t8686) * t8356 - t8439 * t8690));
t8502 = pkin(3) * t8599;
t8517 = t8374 * t8380;
t8521 = t8374 * t8376;
t8565 = t8345 * t8376;
t8045 = (pkin(2) * t8380 - t8345 * t8521) * t8354 + (pkin(2) * t8565 + t8517) * t8351 - t8376 * t8502;
t8561 = t8345 * t8380;
t8046 = (pkin(2) * t8376 + t8345 * t8517) * t8354 + (-pkin(2) * t8561 + t8521) * t8351 + t8380 * t8502;
t7954 = t8045 * t8355 + t8046 * t8352;
t7955 = -t8045 * t8352 + t8046 * t8355;
t8174 = t8351 * t8565 + t8380 * t8354;
t8175 = t8351 * t8561 - t8376 * t8354;
t8067 = t8174 * t8355 - t8175 * t8352;
t8068 = t8174 * t8352 + t8175 * t8355;
t8249 = pkin(1) * t8566;
t8259 = t8317 - pkin(6);
t8133 = t8259 * t8701 + t8249;
t8189 = -pkin(6) * t8344 * t8374 - pkin(1) * t8700;
t8677 = (-((-t8067 * t8330 + t8068 * t8331) * t8311 + (t8067 * t8331 + t8068 * t8330) * t8307) * t8695 + ((t7954 * t8330 + t7955 * t8331) * t8311 - (t7954 * t8331 - t7955 * t8330) * t8307) * t8353 + t7999 * pkin(2) * t8599) / ((pkin(1) * t8502 + t8133 * t8354 + t8189 * t8351) * t8353 + t8426);
t8501 = pkin(3) * t8594;
t8516 = t8374 * t8381;
t8520 = t8374 * t8377;
t8564 = t8345 * t8377;
t8049 = (pkin(2) * t8381 - t8345 * t8520) * t8366 + (pkin(2) * t8564 + t8516) * t8357 - t8377 * t8501;
t8560 = t8345 * t8381;
t8052 = (pkin(2) * t8377 + t8345 * t8516) * t8366 + (-pkin(2) * t8560 + t8520) * t8357 + t8381 * t8501;
t7956 = t8049 * t8367 + t8052 * t8358;
t7959 = -t8049 * t8358 + t8052 * t8367;
t8180 = t8357 * t8564 + t8381 * t8366;
t8183 = t8357 * t8560 - t8377 * t8366;
t8075 = t8180 * t8367 - t8183 * t8358;
t8078 = t8180 * t8358 + t8183 * t8367;
t8273 = t8320 - pkin(6);
t8134 = t8273 * t8701 + t8249;
t8676 = (-((-t8075 * t8330 + t8078 * t8331) * t8312 + (t8075 * t8331 + t8078 * t8330) * t8308) * t8694 + ((t7956 * t8330 + t7959 * t8331) * t8312 - (t7956 * t8331 - t7959 * t8330) * t8308) * t8365 + t8004 * pkin(2) * t8594) / ((pkin(1) * t8501 + t8134 * t8366 + t8189 * t8357) * t8365 + t8425);
t8500 = pkin(3) * t8592;
t8515 = t8374 * t8382;
t8519 = t8374 * t8378;
t8563 = t8345 * t8378;
t8050 = (pkin(2) * t8382 - t8345 * t8519) * t8369 + (pkin(2) * t8563 + t8515) * t8360 - t8378 * t8500;
t8559 = t8345 * t8382;
t8053 = (pkin(2) * t8378 + t8345 * t8515) * t8369 + (-pkin(2) * t8559 + t8519) * t8360 + t8382 * t8500;
t7957 = t8050 * t8370 + t8053 * t8361;
t7960 = -t8050 * t8361 + t8053 * t8370;
t8181 = t8360 * t8563 + t8382 * t8369;
t8184 = t8360 * t8559 - t8378 * t8369;
t8076 = t8181 * t8370 - t8184 * t8361;
t8079 = t8181 * t8361 + t8184 * t8370;
t8275 = t8321 - pkin(6);
t8135 = t8275 * t8701 + t8249;
t8675 = (-((-t8076 * t8330 + t8079 * t8331) * t8313 + (t8076 * t8331 + t8079 * t8330) * t8309) * t8693 + ((t7957 * t8330 + t7960 * t8331) * t8313 - (t7957 * t8331 - t7960 * t8330) * t8309) * t8368 + t8006 * pkin(2) * t8592) / ((pkin(1) * t8500 + t8135 * t8369 + t8189 * t8360) * t8368 + t8424);
t8499 = pkin(3) * t8590;
t8514 = t8374 * t8383;
t8518 = t8374 * t8379;
t8562 = t8345 * t8379;
t8051 = (pkin(2) * t8383 - t8345 * t8518) * t8372 + (pkin(2) * t8562 + t8514) * t8363 - t8379 * t8499;
t8558 = t8345 * t8383;
t8054 = (pkin(2) * t8379 + t8345 * t8514) * t8372 + (-pkin(2) * t8558 + t8518) * t8363 + t8383 * t8499;
t7958 = t8051 * t8373 + t8054 * t8364;
t7961 = -t8051 * t8364 + t8054 * t8373;
t8182 = t8363 * t8562 + t8383 * t8372;
t8185 = t8363 * t8558 - t8379 * t8372;
t8077 = t8182 * t8373 - t8185 * t8364;
t8080 = t8182 * t8364 + t8185 * t8373;
t8277 = t8322 - pkin(6);
t8136 = t8277 * t8701 + t8249;
t8674 = (-((-t8077 * t8330 + t8080 * t8331) * t8314 + (t8077 * t8331 + t8080 * t8330) * t8310) * t8692 + ((t7958 * t8330 + t7961 * t8331) * t8314 - (t7958 * t8331 - t7961 * t8330) * t8310) * t8371 + t8008 * pkin(2) * t8590) / ((pkin(1) * t8499 + t8136 * t8372 + t8189 * t8363) * t8371 + t8423);
t8250 = pkin(1) * t8351 + t8374;
t8329 = pkin(1) * t8374;
t8422 = pkin(1) * t8317 - pkin(6) * t8270;
t8448 = t8345 * t8553;
t8453 = t8344 * t8566;
t8508 = -0.2e1 * (t8345 + 0.1e1) * (t8345 - 0.1e1);
t8384 = pkin(3) ^ 2;
t8605 = t8333 * t8384;
t8609 = t8260 * t8344;
t8613 = (t8318 + t8335) * (t8318 + t8336);
t8673 = (t8329 * t8354 + (-t8345 * t8462 - t8250 + (t8334 * t8508 + t8332) * t8374) * t8260 + ((t8332 * t8613 + t8353 * t8682 - t8604 - t8605) * t8354 - t8258 * t8453) * t8351) / ((t8353 * t8249 - (t8691 - t8699) * t8609) * t8354 - t8448 * t8706 + t8344 * (t8422 * t8353 + (t8270 + pkin(1)) * t8699));
t8254 = pkin(1) * t8357 + t8374;
t8421 = pkin(1) * t8320 - pkin(6) * t8287;
t8447 = t8345 * t8544;
t8603 = t8337 * t8384;
t8608 = t8278 * t8344;
t8612 = (t8323 + t8335) * (t8323 + t8336);
t8672 = (t8329 * t8366 + (-t8345 * t8460 - t8254 + (t8338 * t8508 + t8332) * t8374) * t8278 + ((t8332 * t8612 + t8365 * t8682 - t8603 - t8604) * t8366 - t8272 * t8453) * t8357) / ((t8365 * t8249 - (t8690 - t8698) * t8608) * t8366 - t8447 * t8705 + t8344 * (t8421 * t8365 + (t8287 + pkin(1)) * t8698));
t8255 = pkin(1) * t8360 + t8374;
t8420 = pkin(1) * t8321 - pkin(6) * t8288;
t8446 = t8345 * t8537;
t8602 = t8339 * t8384;
t8607 = t8281 * t8344;
t8611 = (t8325 + t8335) * (t8325 + t8336);
t8671 = (t8329 * t8369 + (-t8345 * t8459 - t8255 + (t8340 * t8508 + t8332) * t8374) * t8281 + ((t8332 * t8611 + t8368 * t8682 - t8602 - t8604) * t8369 - t8274 * t8453) * t8360) / ((t8368 * t8249 - (t8689 - t8697) * t8607) * t8369 - t8446 * t8704 + t8344 * (t8420 * t8368 + (t8288 + pkin(1)) * t8697));
t8256 = pkin(1) * t8363 + t8374;
t8419 = pkin(1) * t8322 - pkin(6) * t8289;
t8445 = t8345 * t8530;
t8601 = t8341 * t8384;
t8606 = t8284 * t8344;
t8610 = (t8327 + t8335) * (t8327 + t8336);
t8670 = (t8329 * t8372 + (-t8345 * t8458 - t8256 + (t8342 * t8508 + t8332) * t8374) * t8284 + ((t8332 * t8610 + t8371 * t8682 - t8601 - t8604) * t8372 - t8276 * t8453) * t8363) / ((t8371 * t8249 - (t8688 - t8696) * t8606) * t8372 - t8445 * t8703 + t8344 * (t8419 * t8371 + (t8289 + pkin(1)) * t8696));
t8452 = t8344 * t8577;
t8669 = ((t8212 * t8600 + t8683) * t8353 + (-pkin(6) * t8452 - t8173 * t8332 + t8198 * t8354) * t8350 + (-(t8301 * t8332 - t8334 + 0.1e1) * t8350 * t8353 + (0.2e1 * t8333 - 0.1e1) * t8452) * pkin(3)) * t8732;
t8451 = t8344 * t8574;
t8668 = ((t8214 * t8600 + t8683) * t8365 + (-pkin(6) * t8451 - t8177 * t8332 + t8199 * t8366) * t8356 + (-(t8302 * t8332 - t8338 + 0.1e1) * t8356 * t8365 + (0.2e1 * t8337 - 0.1e1) * t8451) * pkin(3)) * t8733;
t8450 = t8344 * t8572;
t8667 = ((t8215 * t8600 + t8683) * t8368 + (-pkin(6) * t8450 - t8178 * t8332 + t8200 * t8369) * t8359 + (-(t8303 * t8332 - t8340 + 0.1e1) * t8359 * t8368 + (0.2e1 * t8339 - 0.1e1) * t8450) * pkin(3)) * t8734;
t8449 = t8344 * t8570;
t8666 = ((t8216 * t8600 + t8683) * t8371 + (-pkin(6) * t8449 - t8179 * t8332 + t8201 * t8372) * t8362 + (-(t8304 * t8332 - t8342 + 0.1e1) * t8362 * t8371 + (0.2e1 * t8341 - 0.1e1) * t8449) * pkin(3)) * t8735;
t8137 = t8448 - t8599;
t8549 = t8355 * t8354;
t8109 = -t8137 * t8352 + t8353 * t8549;
t8190 = g(1) * t8307 - t8311 * g(2);
t8194 = g(1) * t8311 + g(2) * t8307;
t8578 = t8345 * t8350;
t7941 = t8109 * t8194 - (t8137 * t8355 + t8352 * t8550) * t8190 + g(3) * (t8344 * t8553 + t8578);
t8386 = pkin(2) ^ 2;
t8435 = t8354 * pkin(6) * t8695;
t8509 = pkin(3) * t8702;
t7949 = 0.1e1 / (t8133 * t8550 + (t8189 * t8353 - t8333 * t8509) * t8351 + (-t8435 + (pkin(2) * t8270 + t8386 * t8354 + t8706) * t8350) * t8344);
t8665 = t7941 * t7949;
t8554 = t8351 * t8352;
t8140 = t8345 * t8554 - t8549;
t8597 = t8344 * t8353;
t8406 = t8140 * t8350 + t8352 * t8597;
t8551 = t8352 * t8354;
t8557 = t8350 * t8351;
t7942 = -t8406 * t8194 - t8190 * ((t8345 * t8557 + t8597) * t8355 + t8350 * t8551) - g(3) * (-t8344 * t8557 + t8576);
t8664 = t7942 * t7949;
t8150 = t8447 - t8594;
t8526 = t8367 * t8366;
t8110 = -t8150 * t8358 + t8365 * t8526;
t8191 = g(1) * t8308 - t8312 * g(2);
t8195 = g(1) * t8312 + g(2) * t8308;
t8575 = t8345 * t8356;
t7943 = t8195 * t8110 - (t8150 * t8367 + t8358 * t8527) * t8191 + g(3) * (t8344 * t8544 + t8575);
t8434 = t8366 * pkin(6) * t8694;
t7950 = 0.1e1 / (t8134 * t8527 + (t8189 * t8365 - t8337 * t8509) * t8357 + (-t8434 + (pkin(2) * t8287 + t8386 * t8366 + t8705) * t8356) * t8344);
t8663 = t7943 * t7950;
t8545 = t8357 * t8358;
t8153 = t8345 * t8545 - t8526;
t8588 = t8344 * t8365;
t8405 = t8153 * t8356 + t8358 * t8588;
t8542 = t8358 * t8366;
t8548 = t8356 * t8357;
t7944 = -t8195 * t8405 - ((t8345 * t8548 + t8588) * t8367 + t8356 * t8542) * t8191 - g(3) * (-t8344 * t8548 + t8569);
t8662 = t7944 * t7950;
t8151 = t8446 - t8592;
t8524 = t8370 * t8369;
t8111 = -t8151 * t8361 + t8368 * t8524;
t8192 = g(1) * t8309 - t8313 * g(2);
t8196 = g(1) * t8313 + g(2) * t8309;
t8573 = t8345 * t8359;
t7945 = t8196 * t8111 - t8192 * (t8151 * t8370 + t8361 * t8525) + g(3) * (t8344 * t8537 + t8573);
t8433 = t8369 * pkin(6) * t8693;
t7951 = 0.1e1 / (t8135 * t8525 + (t8189 * t8368 - t8339 * t8509) * t8360 + (-t8433 + (pkin(2) * t8288 + t8386 * t8369 + t8704) * t8359) * t8344);
t8661 = t7945 * t7951;
t8538 = t8360 * t8361;
t8154 = t8345 * t8538 - t8524;
t8585 = t8344 * t8368;
t8404 = t8154 * t8359 + t8361 * t8585;
t8535 = t8361 * t8369;
t8541 = t8359 * t8360;
t7946 = -t8196 * t8404 - t8192 * ((t8345 * t8541 + t8585) * t8370 + t8359 * t8535) - g(3) * (-t8344 * t8541 + t8568);
t8660 = t7946 * t7951;
t8152 = t8445 - t8590;
t8522 = t8373 * t8372;
t8112 = -t8152 * t8364 + t8371 * t8522;
t8193 = g(1) * t8310 - t8314 * g(2);
t8197 = g(1) * t8314 + g(2) * t8310;
t8571 = t8345 * t8362;
t7947 = t8112 * t8197 - t8193 * (t8152 * t8373 + t8364 * t8523) + g(3) * (t8344 * t8530 + t8571);
t8432 = t8372 * pkin(6) * t8692;
t7952 = 0.1e1 / (t8136 * t8523 + (t8189 * t8371 - t8341 * t8509) * t8363 + (-t8432 + (pkin(2) * t8289 + t8386 * t8372 + t8703) * t8362) * t8344);
t8659 = t7947 * t7952;
t8531 = t8363 * t8364;
t8155 = t8345 * t8531 - t8522;
t8582 = t8344 * t8371;
t8403 = t8155 * t8362 + t8364 * t8582;
t8528 = t8364 * t8372;
t8534 = t8362 * t8363;
t7948 = -t8403 * t8197 - t8193 * ((t8345 * t8534 + t8582) * t8373 + t8362 * t8528) - g(3) * (-t8344 * t8534 + t8567);
t8658 = t7948 * t7952;
t7967 = 0.1e1 / (((t8259 * t8319 + t8422) * t8344 - t8165) * t8353 + t8426);
t8552 = t8351 * t8355;
t8142 = t8345 * t8552 + t8551;
t8595 = t8344 * t8355;
t7968 = -t8194 * (t8142 * t8350 + t8353 * t8595) + t8190 * t8406;
t8657 = t7967 * t7968;
t7969 = t8194 * (t8142 * t8353 - t8350 * t8595) + t8109 * t8190;
t8656 = t7967 * t7969;
t8141 = t8345 * t8549 - t8554;
t8143 = t8345 * t8551 + t8552;
t8019 = t8141 * t8194 - t8143 * t8190;
t8655 = t7967 * t8019;
t8020 = -t8140 * t8190 + t8142 * t8194;
t8654 = t7967 * t8020;
t8033 = -t8344 * t8497 + (-pkin(3) * t8578 - t8212 * t8344) * t8353 - pkin(2) * t8578;
t8653 = t7967 * t8033;
t8086 = t8190 * t8355 + t8194 * t8352;
t8652 = t7967 * t8086;
t8087 = -t8190 * t8352 + t8194 * t8355;
t8651 = t7967 * t8087;
t7974 = 0.1e1 / (((t8277 * t8328 + t8419) * t8344 - t8169) * t8371 + t8423);
t8529 = t8363 * t8373;
t8161 = t8345 * t8529 + t8528;
t8580 = t8344 * t8373;
t7983 = -t8197 * (t8161 * t8362 + t8371 * t8580) + t8193 * t8403;
t8650 = t7974 * t7983;
t7984 = (t8161 * t8371 - t8362 * t8580) * t8197 + t8112 * t8193;
t8649 = t7974 * t7984;
t8160 = t8345 * t8522 - t8531;
t8164 = t8345 * t8528 + t8529;
t8031 = t8160 * t8197 - t8164 * t8193;
t8648 = t7974 * t8031;
t8032 = -t8155 * t8193 + t8161 * t8197;
t8647 = t7974 * t8032;
t8038 = -t8344 * t8491 + (-pkin(3) * t8571 - t8216 * t8344) * t8371 - pkin(2) * t8571;
t8646 = t7974 * t8038;
t8091 = t8193 * t8373 + t8197 * t8364;
t8645 = t7974 * t8091;
t8094 = -t8193 * t8364 + t8197 * t8373;
t8644 = t7974 * t8094;
t7976 = 0.1e1 / (((t8273 * t8324 + t8421) * t8344 - t8167) * t8365 + t8425);
t8543 = t8357 * t8367;
t8157 = t8345 * t8543 + t8542;
t8586 = t8344 * t8367;
t7979 = -t8195 * (t8157 * t8356 + t8365 * t8586) + t8191 * t8405;
t8643 = t7976 * t7979;
t7980 = t8195 * (t8157 * t8365 - t8356 * t8586) + t8110 * t8191;
t8642 = t7976 * t7980;
t8156 = t8345 * t8526 - t8545;
t8162 = t8345 * t8542 + t8543;
t8027 = t8156 * t8195 - t8162 * t8191;
t8641 = t7976 * t8027;
t8028 = -t8153 * t8191 + t8157 * t8195;
t8640 = t7976 * t8028;
t8036 = -t8344 * t8495 + (-pkin(3) * t8575 - t8214 * t8344) * t8365 - pkin(2) * t8575;
t8639 = t7976 * t8036;
t8089 = t8191 * t8367 + t8195 * t8358;
t8638 = t7976 * t8089;
t8092 = -t8191 * t8358 + t8195 * t8367;
t8637 = t7976 * t8092;
t7978 = 0.1e1 / (((t8275 * t8326 + t8420) * t8344 - t8168) * t8368 + t8424);
t8536 = t8360 * t8370;
t8159 = t8345 * t8536 + t8535;
t8583 = t8344 * t8370;
t7981 = -t8196 * (t8159 * t8359 + t8368 * t8583) + t8192 * t8404;
t8636 = t7978 * t7981;
t7982 = t8196 * (t8159 * t8368 - t8359 * t8583) + t8192 * t8111;
t8635 = t7978 * t7982;
t8158 = t8345 * t8524 - t8538;
t8163 = t8345 * t8535 + t8536;
t8029 = t8158 * t8196 - t8163 * t8192;
t8634 = t7978 * t8029;
t8030 = -t8154 * t8192 + t8159 * t8196;
t8633 = t7978 * t8030;
t8037 = -t8344 * t8493 + (-pkin(3) * t8573 - t8215 * t8344) * t8368 - pkin(2) * t8573;
t8632 = t7978 * t8037;
t8090 = t8192 * t8370 + t8196 * t8361;
t8631 = t7978 * t8090;
t8093 = -t8192 * t8361 + t8196 * t8370;
t8630 = t7978 * t8093;
t7986 = 0.1e1 / (t8113 + t8344 * (t8353 * t8410 - t8435 + t8498));
t8225 = g(1) * t8352 - g(2) * t8355;
t8226 = g(1) * t8355 + g(2) * t8352;
t8305 = t8344 * g(3);
t7997 = (t8305 + (-t8225 * t8311 - t8226 * t8307) * t8345) * t8354 + (t8225 * t8307 - t8226 * t8311) * t8351;
t8629 = t7986 * t7997;
t8009 = -g(3) * t8596 + t8141 * t8190 + t8143 * t8194;
t8628 = t7986 * t8009;
t8598 = t8344 * t8351;
t8010 = g(3) * t8598 - t8140 * t8194 - t8142 * t8190;
t8627 = t7986 * t8010;
t7991 = 0.1e1 / (t8114 + t8344 * (t8365 * t8409 - t8434 + t8496));
t8230 = g(1) * t8358 - g(2) * t8367;
t8231 = g(1) * t8367 + g(2) * t8358;
t8000 = (t8305 + (-t8230 * t8312 - t8231 * t8308) * t8345) * t8366 + t8357 * (t8230 * t8308 - t8231 * t8312);
t8626 = t7991 * t8000;
t8011 = -g(3) * t8587 + t8156 * t8191 + t8162 * t8195;
t8625 = t7991 * t8011;
t8593 = t8344 * t8357;
t8012 = g(3) * t8593 - t8195 * t8153 - t8157 * t8191;
t8624 = t7991 * t8012;
t7992 = 0.1e1 / (t8115 + t8344 * (t8368 * t8408 - t8433 + t8494));
t8232 = g(1) * t8361 - g(2) * t8370;
t8233 = g(1) * t8370 + g(2) * t8361;
t8001 = (t8305 + (-t8232 * t8313 - t8233 * t8309) * t8345) * t8369 + t8360 * (t8232 * t8309 - t8233 * t8313);
t8623 = t7992 * t8001;
t8013 = -g(3) * t8584 + t8158 * t8192 + t8163 * t8196;
t8622 = t7992 * t8013;
t8591 = t8344 * t8360;
t8014 = g(3) * t8591 - t8196 * t8154 - t8159 * t8192;
t8621 = t7992 * t8014;
t7993 = 0.1e1 / (t8116 + t8344 * (t8371 * t8407 - t8432 + t8492));
t8234 = g(1) * t8364 - g(2) * t8373;
t8235 = g(1) * t8373 + g(2) * t8364;
t8002 = (t8305 + (-t8234 * t8314 - t8235 * t8310) * t8345) * t8372 + (t8234 * t8310 - t8235 * t8314) * t8363;
t8620 = t7993 * t8002;
t8015 = -g(3) * t8581 + t8160 * t8193 + t8164 * t8197;
t8619 = t7993 * t8015;
t8589 = t8344 * t8363;
t8016 = g(3) * t8589 - t8155 * t8197 - t8161 * t8193;
t8618 = t7993 * t8016;
t8617 = (pkin(1) + 0.2e1 * t8270) * t8260;
t8616 = (pkin(1) + 0.2e1 * t8287) * t8278;
t8615 = (pkin(1) + 0.2e1 * t8288) * t8281;
t8614 = (pkin(1) + 0.2e1 * t8289) * t8284;
t8138 = t8307 * t8355 + t8311 * t8352;
t8139 = -t8307 * t8352 + t8311 * t8355;
t8556 = t8351 * (t8138 * t8700 - t8139 * t8374);
t8555 = t8351 * (t8138 * t8374 + t8139 * t8700);
t8144 = t8308 * t8367 + t8312 * t8358;
t8145 = -t8308 * t8358 + t8312 * t8367;
t8547 = t8357 * (t8144 * t8700 - t8145 * t8374);
t8546 = t8357 * (t8144 * t8374 + t8145 * t8700);
t8146 = t8309 * t8370 + t8313 * t8361;
t8147 = -t8309 * t8361 + t8313 * t8370;
t8540 = t8360 * (t8146 * t8700 - t8147 * t8374);
t8539 = t8360 * (t8146 * t8374 + t8147 * t8700);
t8148 = t8310 * t8373 + t8314 * t8364;
t8149 = -t8310 * t8364 + t8314 * t8373;
t8533 = t8363 * (t8148 * t8700 - t8149 * t8374);
t8532 = t8363 * (t8148 * t8374 + t8149 * t8700);
t8486 = t7997 * t8681;
t8485 = t8001 * t8680;
t8484 = t8002 * t8679;
t8483 = t8000 * t8678;
t8482 = t7997 * t8669;
t8481 = t8000 * t8668;
t8480 = t8001 * t8667;
t8479 = t8002 * t8666;
t8478 = t8350 * t8629;
t8477 = t8353 * t8629;
t8476 = t8356 * t8626;
t8475 = t8365 * t8626;
t8474 = t8359 * t8623;
t8473 = t8368 * t8623;
t8472 = t8362 * t8620;
t8471 = t8371 * t8620;
t8470 = t8138 * t8599;
t8469 = t8139 * t8599;
t8468 = t8144 * t8594;
t8467 = t8145 * t8594;
t8466 = t8146 * t8592;
t8465 = t8147 * t8592;
t8464 = t8148 * t8590;
t8463 = t8149 * t8590;
t8461 = t8260 * t8566;
t8457 = t8278 * t8566;
t8456 = t8281 * t8566;
t8455 = t8284 * t8566;
t8444 = t8352 * t8613;
t8443 = t8358 * t8612;
t8442 = t8361 * t8611;
t8441 = t8364 * t8610;
t8430 = t8355 * t8461;
t8429 = t8367 * t8457;
t8428 = t8370 * t8456;
t8427 = t8373 * t8455;
t8261 = 0.2e1 * t8319 + pkin(1);
t8343 = t8374 ^ 2;
t8418 = pkin(1) * t8319 + t8261 * t8270 + t8334 * t8604 + t8343;
t8279 = 0.2e1 * t8324 + pkin(1);
t8417 = pkin(1) * t8324 + t8279 * t8287 + t8338 * t8604 + t8343;
t8282 = 0.2e1 * t8326 + pkin(1);
t8416 = pkin(1) * t8326 + t8282 * t8288 + t8340 * t8604 + t8343;
t8285 = 0.2e1 * t8328 + pkin(1);
t8415 = pkin(1) * t8328 + t8285 * t8289 + t8342 * t8604 + t8343;
t8394 = ((t8301 * t8317 - pkin(6)) * t8353 + t8350 * t8173) * t8600 - (t8351 * t8490 + t8513) * t8332 - t8351 * t8716;
t8393 = ((t8302 * t8320 - pkin(6)) * t8365 + t8356 * t8177) * t8600 - (t8357 * t8489 + t8512) * t8332 - t8357 * t8717;
t8392 = ((t8303 * t8321 - pkin(6)) * t8368 + t8359 * t8178) * t8600 - (t8360 * t8488 + t8511) * t8332 - t8360 * t8718;
t8391 = ((t8304 * t8322 - pkin(6)) * t8371 + t8362 * t8179) * t8600 - (t8363 * t8487 + t8510) * t8332 - t8363 * t8719;
t8203 = g(1) * t8331 + g(2) * t8330;
t8202 = g(1) * t8330 - g(2) * t8331;
t8120 = -t8276 * t8589 + t8284 * t8345;
t8119 = -t8274 * t8591 + t8281 * t8345;
t8118 = -t8272 * t8593 + t8278 * t8345;
t8117 = -t8258 * t8598 + t8260 * t8345;
t8084 = 0.2e1 * t8364 * t8455 + t8373 * t8610;
t8083 = 0.2e1 * t8361 * t8456 + t8370 * t8611;
t8082 = 0.2e1 * t8358 * t8457 + t8367 * t8612;
t8081 = 0.2e1 * t8352 * t8461 + t8355 * t8613;
t8074 = t8120 * t8373 + t8256 * t8364;
t8073 = t8119 * t8370 + t8255 * t8361;
t8072 = t8118 * t8367 + t8254 * t8358;
t8071 = -t8120 * t8364 + t8256 * t8373;
t8070 = -t8119 * t8361 + t8255 * t8370;
t8069 = -t8118 * t8358 + t8254 * t8367;
t8066 = t8117 * t8355 + t8250 * t8352;
t8065 = -t8117 * t8352 + t8250 * t8355;
t8064 = -t8276 * t8606 + t8570 * t8610;
t8063 = -t8274 * t8607 + t8572 * t8611;
t8062 = -t8272 * t8608 + t8574 * t8612;
t8061 = -t8258 * t8609 + t8577 * t8613;
t8044 = -t8148 * t8570 + t8149 * t8372;
t8043 = t8148 * t8372 + t8149 * t8570;
t8042 = -t8146 * t8572 + t8147 * t8369;
t8041 = t8146 * t8369 + t8147 * t8572;
t8040 = -t8144 * t8574 + t8145 * t8366;
t8039 = t8144 * t8366 + t8145 * t8574;
t8035 = -t8138 * t8577 + t8139 * t8354;
t8034 = t8138 * t8354 + t8139 * t8577;
t8026 = t8064 * t8373 + t8364 * t8614;
t8025 = t8063 * t8370 + t8361 * t8615;
t8024 = t8062 * t8367 + t8358 * t8616;
t8023 = -t8064 * t8364 + t8373 * t8614;
t8022 = -t8063 * t8361 + t8370 * t8615;
t8021 = -t8062 * t8358 + t8367 * t8616;
t8018 = t8061 * t8355 + t8352 * t8617;
t8017 = -t8061 * t8352 + t8355 * t8617;
t7932 = -t8044 * t8692 + (-pkin(3) * t8464 + (-pkin(2) * t8149 - t8148 * t8566) * t8372 + t8533) * t8371 - pkin(2) * t8464;
t7931 = -t8043 * t8692 + (pkin(3) * t8463 + (-pkin(2) * t8148 + t8149 * t8566) * t8372 - t8532) * t8371 + pkin(2) * t8463;
t7930 = -t8042 * t8693 + (-pkin(3) * t8466 + (-pkin(2) * t8147 - t8146 * t8566) * t8369 + t8540) * t8368 - pkin(2) * t8466;
t7929 = -t8041 * t8693 + (pkin(3) * t8465 + (-pkin(2) * t8146 + t8147 * t8566) * t8369 - t8539) * t8368 + pkin(2) * t8465;
t7928 = -t8040 * t8694 + (-pkin(3) * t8468 + (-pkin(2) * t8145 - t8144 * t8566) * t8366 + t8547) * t8365 - pkin(2) * t8468;
t7927 = -t8039 * t8694 + (pkin(3) * t8467 + (-pkin(2) * t8144 + t8145 * t8566) * t8366 - t8546) * t8365 + pkin(2) * t8467;
t7926 = -t8035 * t8695 + (-pkin(3) * t8470 + (-pkin(2) * t8139 - t8138 * t8566) * t8354 + t8556) * t8353 - pkin(2) * t8470;
t7925 = -t8034 * t8695 + (pkin(3) * t8469 + (-pkin(2) * t8138 + t8139 * t8566) * t8354 - t8555) * t8353 + pkin(2) * t8469;
t7924 = t8148 * t8399 + t8149 * t8391;
t7923 = -t8148 * t8391 + t8149 * t8399;
t7922 = t8146 * t8400 + t8147 * t8392;
t7921 = -t8146 * t8392 + t8147 * t8400;
t7920 = t8144 * t8401 + t8145 * t8393;
t7919 = -t8144 * t8393 + t8145 * t8401;
t7918 = t8138 * t8402 + t8139 * t8394;
t7917 = -t8138 * t8394 + t8139 * t8402;
t7916 = ((-0.2e1 * t8427 + t8441) * t8314 + t8084 * t8310) * t8342 + (t8023 * t8310 + t8026 * t8314) * t8372 + t8374 * (t8071 * t8310 + t8074 * t8314);
t7915 = ((-0.2e1 * t8428 + t8442) * t8313 + t8083 * t8309) * t8340 + (t8022 * t8309 + t8025 * t8313) * t8369 + t8374 * (t8070 * t8309 + t8073 * t8313);
t7914 = ((-0.2e1 * t8429 + t8443) * t8312 + t8082 * t8308) * t8338 + (t8021 * t8308 + t8024 * t8312) * t8366 + t8374 * (t8069 * t8308 + t8072 * t8312);
t7913 = (t8084 * t8314 + 0.2e1 * (t8427 - t8441 / 0.2e1) * t8310) * t8342 + (t8023 * t8314 - t8026 * t8310) * t8372 + (t8071 * t8314 - t8074 * t8310) * t8374;
t7912 = (t8083 * t8313 + 0.2e1 * (t8428 - t8442 / 0.2e1) * t8309) * t8340 + (t8022 * t8313 - t8025 * t8309) * t8369 + (t8070 * t8313 - t8073 * t8309) * t8374;
t7911 = (t8082 * t8312 + 0.2e1 * (t8429 - t8443 / 0.2e1) * t8308) * t8338 + (t8021 * t8312 - t8024 * t8308) * t8366 + (t8069 * t8312 - t8072 * t8308) * t8374;
t7910 = ((-0.2e1 * t8430 + t8444) * t8311 + t8081 * t8307) * t8334 + (t8017 * t8307 + t8018 * t8311) * t8354 + t8374 * (t8065 * t8307 + t8066 * t8311);
t7909 = (t8081 * t8311 + 0.2e1 * (t8430 - t8444 / 0.2e1) * t8307) * t8334 + (t8017 * t8311 - t8018 * t8307) * t8354 + (t8065 * t8311 - t8066 * t8307) * t8374;
t7900 = (-(t8330 * t8383 + t8331 * t8379) * (t8415 * t8149 + (t8044 * t8601 + (t8149 * t8285 - 0.2e1 * t8533) * t8327) * t8372 + t8708 * t8148) + (t8415 * t8148 + (t8043 * t8601 + (t8148 * t8285 + 0.2e1 * t8532) * t8327) * t8372 - t8708 * t8149) * (-t8330 * t8379 + t8331 * t8383)) * t8735 * t8579;
t7899 = (-(t8330 * t8382 + t8331 * t8378) * (t8416 * t8147 + (t8042 * t8602 + (t8147 * t8282 - 0.2e1 * t8540) * t8325) * t8369 + t8709 * t8146) + (t8416 * t8146 + (t8041 * t8602 + (t8146 * t8282 + 0.2e1 * t8539) * t8325) * t8369 - t8709 * t8147) * (-t8330 * t8378 + t8331 * t8382)) * t8734 * t8579;
t7898 = (-(t8417 * t8145 + (t8040 * t8603 + (t8145 * t8279 - 0.2e1 * t8547) * t8323) * t8366 + t8710 * t8144) * (t8330 * t8381 + t8331 * t8377) + (t8417 * t8144 + (t8039 * t8603 + (t8144 * t8279 + 0.2e1 * t8546) * t8323) * t8366 - t8710 * t8145) * (-t8330 * t8377 + t8331 * t8381)) * t8733 * t8579;
t7897 = (-(t8330 * t8380 + t8331 * t8376) * (t8418 * t8139 + (t8035 * t8605 + (t8139 * t8261 - 0.2e1 * t8556) * t8318) * t8354 + t8711 * t8138) + (-t8330 * t8376 + t8331 * t8380) * (t8418 * t8138 + (t8034 * t8605 + (t8138 * t8261 + 0.2e1 * t8555) * t8318) * t8354 - t8711 * t8139)) * t8732 * t8579;
t1 = [0, t7926 * t8652 + t7928 * t8638 + t7930 * t8631 + t7932 * t8645, t7926 * t8651 + t7928 * t8637 + t7930 * t8630 + t7932 * t8644, 0, 0, 0, 0, 0, t7917 * t8628 + t7919 * t8625 + t7921 * t8622 + t7923 * t8619 + t7926 * t8654 + t7928 * t8640 + t7930 * t8633 + t7932 * t8647, t7917 * t8627 + t7919 * t8624 + t7921 * t8621 + t7923 * t8618 + t7926 * t8655 + t7928 * t8641 + t7930 * t8634 + t7932 * t8648, 0, 0, 0, 0, 0, -t7917 * t8477 - t7919 * t8475 - t7921 * t8473 - t7923 * t8471 + t7926 * t8656 + t7928 * t8642 + t7930 * t8635 + t7932 * t8649 + (-t7909 * t8664 - t7911 * t8662 - t7912 * t8660 - t7913 * t8658) * t8579, t7917 * t8478 + t7919 * t8476 + t7921 * t8474 + t7923 * t8472 + t7926 * t8657 + t7928 * t8643 + t7930 * t8636 + t7932 * t8650 + (-t7909 * t8665 - t7911 * t8663 - t7912 * t8661 - t7913 * t8659) * t8579, 0, 0, 0, -t8202 * t8330 - t8203 * t8331; 0, t7925 * t8652 + t7927 * t8638 + t7929 * t8631 + t7931 * t8645, t7925 * t8651 + t7927 * t8637 + t7929 * t8630 + t7931 * t8644, 0, 0, 0, 0, 0, t7918 * t8628 + t7920 * t8625 + t7922 * t8622 + t7924 * t8619 + t7925 * t8654 + t7927 * t8640 + t7929 * t8633 + t7931 * t8647, t7918 * t8627 + t7920 * t8624 + t7922 * t8621 + t7924 * t8618 + t7925 * t8655 + t7927 * t8641 + t7929 * t8634 + t7931 * t8648, 0, 0, 0, 0, 0, -t7918 * t8477 - t7920 * t8475 - t7922 * t8473 - t7924 * t8471 + t7925 * t8656 + t7927 * t8642 + t7929 * t8635 + t7931 * t8649 + (-t7910 * t8664 - t7914 * t8662 - t7915 * t8660 - t7916 * t8658) * t8579, t7918 * t8478 + t7920 * t8476 + t7922 * t8474 + t7924 * t8472 + t7925 * t8657 + t7927 * t8643 + t7929 * t8636 + t7931 * t8650 + (-t7910 * t8665 - t7914 * t8663 - t7915 * t8661 - t7916 * t8659) * t8579, 0, 0, 0, t8202 * t8331 - t8203 * t8330; 0, t8033 * t8652 + t8036 * t8638 + t8037 * t8631 + t8038 * t8645, t8033 * t8651 + t8036 * t8637 + t8037 * t8630 + t8038 * t8644, 0, 0, 0, 0, 0, -t8009 * t8669 - t8011 * t8668 - t8013 * t8667 - t8015 * t8666 + t8020 * t8653 + t8028 * t8639 + t8030 * t8632 + t8032 * t8646, -t8010 * t8669 - t8012 * t8668 - t8014 * t8667 - t8016 * t8666 + t8019 * t8653 + t8027 * t8639 + t8029 * t8632 + t8031 * t8646, 0, 0, 0, 0, 0, t8353 * t8482 + t8365 * t8481 + t8368 * t8480 + t8371 * t8479 + t7969 * t8653 + t7984 * t8646 + t7980 * t8639 + t7982 * t8632 + (t7942 * t8673 + t7944 * t8672 + t7946 * t8671 + t7948 * t8670) * t8385, -t8350 * t8482 - t8356 * t8481 - t8359 * t8480 - t8362 * t8479 + t7968 * t8653 + t7983 * t8646 + t7979 * t8639 + t7981 * t8632 + (t7941 * t8673 + t7943 * t8672 + t7945 * t8671 + t7947 * t8670) * t8385, 0, 0, 0, -g(3); 0, t8086 * t8677 + t8089 * t8676 + t8090 * t8675 + t8091 * t8674, t8087 * t8677 + t8092 * t8676 + t8093 * t8675 + t8094 * t8674, 0, 0, 0, 0, 0, t8009 * t8681 + t8011 * t8678 + t8013 * t8680 + t8015 * t8679 + t8020 * t8677 + t8028 * t8676 + t8030 * t8675 + t8032 * t8674, t8010 * t8681 + t8012 * t8678 + t8014 * t8680 + t8016 * t8679 + t8019 * t8677 + t8027 * t8676 + t8029 * t8675 + t8031 * t8674, 0, 0, 0, 0, 0, t7897 * t7942 + t7898 * t7944 + t7899 * t7946 + t7900 * t7948 + t7969 * t8677 + t7980 * t8676 + t7982 * t8675 + t7984 * t8674 - t8353 * t8486 - t8365 * t8483 - t8368 * t8485 - t8371 * t8484, t7897 * t7941 + t7898 * t7943 + t7899 * t7945 + t7900 * t7947 + t7968 * t8677 + t7979 * t8676 + t7981 * t8675 + t7983 * t8674 + t8350 * t8486 + t8356 * t8483 + t8359 * t8485 + t8362 * t8484, 0, t8202, t8203, 0;];
tau_reg  = t1;
