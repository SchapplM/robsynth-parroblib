% Calculate minimal parameter regressor of Gravitation load for parallel robot
% P3PRRR2G3A0
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
% pkin [2x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a3,a4]';
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
% Datum: 2020-03-09 21:20
% Revision: 0f11fd83bca0a8cdff505979e09e2c4d81033460 (2020-02-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = P3PRRR2G3A0_gravload_para_pf_regmin(xP, qJ, g, legFrame, ...
  koppelP, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,3),zeros(3,1),zeros(3,3),zeros(3,3),zeros(2,1)}
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3PRRR2G3A0_gravload_para_pf_regmin: xP has to be [3x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 3]), ...
  'P3PRRR2G3A0_gravload_para_pf_regmin: qJ has to be [3x3] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'P3PRRR2G3A0_gravload_para_pf_regmin: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [2 1]), ...
  'P3PRRR2G3A0_gravload_para_pf_regmin: pkin has to be [2x1] (double)');
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3PRRR2G3A0_gravload_para_pf_regmin: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3PRRR2G3A0_gravload_para_pf_regmin: Koppelpunkt has to be [3x3] (double)');

%% Symbolic Calculation
% From invdyn_para_plfcoord_taugreg_matlab.m
% OptimizationMode: 2
% StartTime: 2020-03-09 21:20:07
% EndTime: 2020-03-09 21:20:07
% DurationCPUTime: 0.16s
% Computational Cost: add. (283->64), mult. (277->109), div. (84->5), fcn. (300->33), ass. (0->72)
t500 = legFrame(3,2);
t525 = sin(t500);
t501 = legFrame(2,2);
t524 = sin(t501);
t502 = legFrame(1,2);
t523 = sin(t502);
t482 = -t500 + qJ(2,3);
t476 = qJ(3,3) + t482;
t470 = sin(t476);
t494 = 0.1e1 / sin(qJ(3,3));
t522 = (pkin(2) * t470 + pkin(1) * sin(t482)) * t494;
t483 = -t501 + qJ(2,2);
t477 = qJ(3,2) + t483;
t471 = sin(t477);
t495 = 0.1e1 / sin(qJ(3,2));
t521 = (pkin(2) * t471 + pkin(1) * sin(t483)) * t495;
t484 = -t502 + qJ(2,1);
t478 = qJ(3,1) + t484;
t472 = sin(t478);
t496 = 0.1e1 / sin(qJ(3,1));
t520 = (pkin(2) * t472 + pkin(1) * sin(t484)) * t496;
t473 = cos(t476);
t519 = (-pkin(2) * t473 - pkin(1) * cos(t482)) * t494;
t474 = cos(t477);
t518 = (-pkin(2) * t474 - pkin(1) * cos(t483)) * t495;
t475 = cos(t478);
t517 = (-pkin(2) * t475 - pkin(1) * cos(t484)) * t496;
t516 = t470 * t494;
t515 = t471 * t495;
t514 = t472 * t496;
t513 = t473 * t494;
t512 = t474 * t495;
t511 = t475 * t496;
t510 = 0.1e1 / pkin(1);
t509 = 0.1e1 / pkin(2);
t508 = cos(qJ(2,1));
t507 = cos(qJ(2,2));
t506 = cos(qJ(2,3));
t505 = sin(qJ(2,1));
t504 = sin(qJ(2,2));
t503 = sin(qJ(2,3));
t499 = qJ(2,1) + qJ(3,1);
t498 = qJ(2,2) + qJ(3,2);
t497 = qJ(2,3) + qJ(3,3);
t493 = cos(t502);
t492 = cos(t501);
t491 = cos(t500);
t490 = cos(t499);
t489 = cos(t498);
t488 = cos(t497);
t487 = sin(t499);
t486 = sin(t498);
t485 = sin(t497);
t469 = t493 * g(1) - t523 * g(2);
t468 = t492 * g(1) - t524 * g(2);
t467 = t491 * g(1) - t525 * g(2);
t466 = t523 * g(1) + t493 * g(2);
t465 = t524 * g(1) + t492 * g(2);
t464 = t525 * g(1) + t491 * g(2);
t457 = t466 * t508 - t469 * t505;
t456 = t466 * t505 + t469 * t508;
t455 = t465 * t507 - t468 * t504;
t454 = t465 * t504 + t468 * t507;
t453 = t464 * t506 - t467 * t503;
t452 = t464 * t503 + t467 * t506;
t451 = t466 * t490 - t469 * t487;
t450 = t466 * t487 + t469 * t490;
t449 = t465 * t489 - t468 * t486;
t448 = t465 * t486 + t468 * t489;
t447 = t464 * t488 - t467 * t485;
t446 = t464 * t485 + t467 * t488;
t1 = [0, 0, (-t452 * t516 - t454 * t515 - t456 * t514) * t510, (-t453 * t516 - t455 * t515 - t457 * t514) * t510, 0, (-t446 * t516 - t448 * t515 - t450 * t514 + (t446 * t522 + t448 * t521 + t450 * t520) * t509) * t510, (-t447 * t516 - t449 * t515 - t451 * t514 + (t447 * t522 + t449 * t521 + t451 * t520) * t509) * t510, -g(1); 0, 0, (t452 * t513 + t454 * t512 + t456 * t511) * t510, (t453 * t513 + t455 * t512 + t457 * t511) * t510, 0, (t446 * t513 + t448 * t512 + t450 * t511 + (t446 * t519 + t448 * t518 + t450 * t517) * t509) * t510, (t447 * t513 + t449 * t512 + t451 * t511 + (t447 * t519 + t449 * t518 + t451 * t517) * t509) * t510, -g(2); -3 * g(3), 0, 0, 0, 0, 0, 0, -g(3);];
tau_reg  = t1;
