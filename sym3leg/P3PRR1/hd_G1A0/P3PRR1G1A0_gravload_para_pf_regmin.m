% Calculate minimal parameter regressor of Gravitation load for parallel robot
% P3PRR1G1A0
% Use Code from Maple symbolic Code Generation
%
% Input:
% xP [3x1]
%   Generalized platform coordinates
% qJ [2x3]
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
% pkin [4x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,d2,d3]';
% koppelP [3x3]
%   coordinates of the platform coupling joints
%   (joints that link the end of legs with platform)
%   in platform coordinates
%   rows: number of leg
%   columns: x-, y-, z-coordinates

% Output:
% tau_reg [3x8]
%   minimal parameter regressor of Gravitation load for parallel robot
%   in task space

% Quelle: HybrDyn-Toolbox
% Datum: 2019-05-03 14:47
% Revision: abbb0d669c4fc7889a31e0cf750ab51a4f2eb1ce (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = P3PRR1G1A0_gravload_para_pf_regmin(xP, qJ, g, legFrame, ...
  koppelP, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(2,3),zeros(3,1),zeros(3,3),zeros(3,3),zeros(4,1)}
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3PRR1G1A0_gravload_para_pf_regmin: xP has to be [3x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [2 3]), ...
  'P3PRR1G1A0_gravload_para_pf_regmin: qJ has to be [2x3] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'P3PRR1G1A0_gravload_para_pf_regmin: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [4 1]), ...
  'P3PRR1G1A0_gravload_para_pf_regmin: pkin has to be [4x1] (double)');
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3PRR1G1A0_gravload_para_pf_regmin: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3PRR1G1A0_gravload_para_pf_regmin: Koppelpunkt has to be [3x3] (double)');

%% Symbolic Calculation
% From invdyn_para_plfcoord_taugreg_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-03 14:47:40
% EndTime: 2019-05-03 14:47:40
% DurationCPUTime: 0.26s
% Computational Cost: add. (134->55), mult. (271->119), div. (45->4), fcn. (295->14), ass. (0->64)
t399 = xP(3);
t385 = sin(t399);
t386 = cos(t399);
t400 = koppelP(3,2);
t403 = koppelP(3,1);
t365 = t385 * t403 + t386 * t400;
t390 = legFrame(3,3);
t379 = sin(t390);
t382 = cos(t390);
t409 = t385 * t400 - t386 * t403;
t427 = t365 * t382 + t409 * t379;
t401 = koppelP(2,2);
t404 = koppelP(2,1);
t366 = t385 * t404 + t386 * t401;
t391 = legFrame(2,3);
t380 = sin(t391);
t383 = cos(t391);
t408 = t385 * t401 - t386 * t404;
t426 = t366 * t383 + t408 * t380;
t402 = koppelP(1,2);
t405 = koppelP(1,1);
t367 = t385 * t405 + t386 * t402;
t392 = legFrame(1,3);
t381 = sin(t392);
t384 = cos(t392);
t407 = t385 * t402 - t386 * t405;
t425 = t367 * t384 + t407 * t381;
t395 = sin(qJ(2,1));
t389 = 0.1e1 / t395;
t424 = t425 * t389;
t393 = sin(qJ(2,3));
t387 = 0.1e1 / t393;
t423 = t427 * t387;
t394 = sin(qJ(2,2));
t388 = 0.1e1 / t394;
t422 = t426 * t388;
t371 = t379 * g(1) - t382 * g(2);
t418 = t371 * t387;
t372 = t380 * g(1) - t383 * g(2);
t417 = t372 * t388;
t373 = t381 * g(1) - t384 * g(2);
t416 = t373 * t389;
t415 = t379 * t387;
t414 = t380 * t388;
t413 = t381 * t389;
t412 = t382 * t387;
t411 = t383 * t388;
t410 = t384 * t389;
t406 = 0.1e1 / pkin(2);
t398 = cos(qJ(2,1));
t397 = cos(qJ(2,2));
t396 = cos(qJ(2,3));
t378 = t386 * g(1) + t385 * g(2);
t377 = t385 * g(1) - t386 * g(2);
t376 = t384 * g(1) + t381 * g(2);
t375 = t383 * g(1) + t380 * g(2);
t374 = t382 * g(1) + t379 * g(2);
t364 = -t373 * t395 + t376 * t398;
t363 = -t372 * t394 + t375 * t397;
t362 = -t371 * t393 + t374 * t396;
t361 = t373 * t398 + t376 * t395;
t360 = t372 * t397 + t375 * t394;
t359 = t371 * t396 + t374 * t393;
t1 = [(-t381 * t395 + t384 * t398) * t416 + (-t380 * t394 + t383 * t397) * t417 + (-t379 * t393 + t382 * t396) * t418, 0, (-t359 * t412 - t360 * t411 - t361 * t410) * t406, (-t362 * t412 - t363 * t411 - t364 * t410) * t406, 0, 0, 0, -t385 * t377 - t386 * t378; (t381 * t398 + t395 * t384) * t416 + (t380 * t397 + t394 * t383) * t417 + (t379 * t396 + t393 * t382) * t418, 0, (-t359 * t415 - t360 * t414 - t361 * t413) * t406, (-t362 * t415 - t363 * t414 - t364 * t413) * t406, 0, 0, 0, t386 * t377 - t385 * t378; ((t381 * t367 - t384 * t407) * t395 - t425 * t398) * t416 + ((t380 * t366 - t383 * t408) * t394 - t426 * t397) * t417 + ((t379 * t365 - t382 * t409) * t393 - t427 * t396) * t418, 0, (t359 * t423 + t360 * t422 + t361 * t424) * t406, (t362 * t423 + t363 * t422 + t364 * t424) * t406, 0, t377, t378, 0;];
tau_reg  = t1;
