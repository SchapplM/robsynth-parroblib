% Calculate vector of centrifugal and coriolis load on the joints for
% P3RPRR1G3P3A0
% Use Code from Maple symbolic Code Generation
%
% Input:
% xP [3x1]
%   Generalized platform coordinates
% xDP [3x1]
%   Generalized platform velocities
% qJ [3x3]
%   Generalized joint coordinates (joint angles)
%   rows: links of the robot
%   columns: number of leg
% legFrame [3x3]
%   base frame orientation for each leg
%   row: number of leg
%   column: Euler angles for the orientation.
%   Euler angle convention from robot definition ("leg_frame")
% koppelP [3x3]
%   coordinates of the platform coupling joints
%   (joints that link the end of legs with platform)
%   in platform coordinates
%   rows: number of leg
%   columns: x-, y-, z-coordinates
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d3,d4,theta2]';
% m [4x1]
%   mass of all robot links (leg links until cut joint, platform)
% mrSges [4x3]
%   first moment of all robot links (mass times center of mass in body frames)
%   rows: links of the robot (leg links until cut joint, platform)
%   columns: x-, y-, z-coordinates
% Ifges [4x6]
%   inertia of all robot links about their respective body frame origins, in body frames
%   rows: links of the robot (leg links until cut joint, platform)
%   columns: xx, yy, zz, xy, xz, yz (see inertial_parameters_convert_par1_par2.m)
%
% Output:
% taucX [3x1]
%   forces required to compensate Coriolis and centrifugal joint torques
%   in platform coordinates

% Quelle: HybrDyn-Toolbox
% Datum: 2020-03-09 21:27
% Revision: 0f11fd83bca0a8cdff505979e09e2c4d81033460 (2020-02-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taucX = P3RPRR1G3P3A0_coriolisvec_para_pf_slag_vp2(xP, xDP, qJ, legFrame, ...
  koppelP, pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,1),zeros(3,3),zeros(3,3),zeros(3,3),zeros(7,1),zeros(3+1,1),zeros(3+1,3),zeros(3+1,6)}
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3RPRR1G3P3A0_coriolisvec_para_pf_slag_vp2: xP has to be [3x1] (double)');
assert(isreal(xDP) && all(size(xDP) == [3 1]), ...
  'P3RPRR1G3P3A0_coriolisvec_para_pf_slag_vp2: xDP has to be [3x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 3]), ...
  'P3RPRR1G3P3A0_coriolisvec_para_pf_slag_vp2: qJ has to be [3x3] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'P3RPRR1G3P3A0_coriolisvec_para_pf_slag_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [4 1]), ...
  'P3RPRR1G3P3A0_coriolisvec_para_pf_slag_vp2: m has to be [4x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [4,3]), ...
  'P3RPRR1G3P3A0_coriolisvec_para_pf_slag_vp2: mrSges has to be [4x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [4 6]), ...
  'P3RPRR1G3P3A0_coriolisvec_para_pf_slag_vp2: Ifges has to be [4x6] (double)'); 
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3RPRR1G3P3A0_coriolisvec_para_pf_slag_vp2: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3RPRR1G3P3A0_coriolisvec_para_pf_slag_vp2: Koppelpunkt has to be [3x3] (double)');

%% Symbolic Calculation
% From coriolisvec_para_plfcoord_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2020-03-09 21:26:52
% EndTime: 2020-03-09 21:26:53
% DurationCPUTime: 1.22s
% Computational Cost: add. (14520->199), mult. (9417->263), div. (1818->7), fcn. (6798->59), ass. (0->142)
t555 = 2 * pkin(2);
t527 = 0.1e1 / pkin(3);
t565 = t527 / 0.2e1;
t513 = legFrame(1,2);
t554 = qJ(1,1) + pkin(7);
t490 = t513 + t554;
t484 = qJ(3,1) + t490;
t491 = -t513 + t554;
t485 = qJ(3,1) + t491;
t462 = cos(t485) + cos(t484);
t503 = qJ(1,1) + t513;
t504 = qJ(1,1) - t513;
t438 = -t462 * pkin(3) + (-cos(t491) - cos(t490)) * pkin(2) + (-cos(t504) - cos(t503)) * pkin(1);
t508 = pkin(7) + qJ(3,1);
t516 = sin(qJ(3,1));
t467 = 0.1e1 / (pkin(1) * sin(t508) + t516 * pkin(2));
t523 = xDP(1);
t550 = t523 * t565;
t432 = t438 * t467 * t550;
t459 = -sin(t484) + sin(t485);
t435 = t459 * pkin(3) + (-sin(t490) + sin(t491)) * pkin(2) + (-sin(t503) + sin(t504)) * pkin(1);
t494 = sin(qJ(1,1) + t508);
t456 = pkin(3) * t494 + pkin(2) * sin(t554) + sin(qJ(1,1)) * pkin(1);
t521 = xDP(3);
t564 = t521 * t527;
t453 = t456 * t467 * t564;
t522 = xDP(2);
t551 = t522 * t565;
t543 = -t551 / 0.2e1;
t566 = t523 / 0.2e1;
t567 = t522 / 0.2e1;
t557 = (t459 * t567 + t462 * t566) * t467;
t568 = t494 * t521;
t417 = t432 / 0.2e1 + t453 / 0.2e1 + (t435 * t543 - t568) * t467 + t557;
t544 = t435 * t551;
t560 = t432 + t453;
t420 = (-t544 - t568) * t467 + t557 + t560;
t429 = -t467 * t568 + t557;
t510 = cos(pkin(7));
t495 = t510 * pkin(1) + pkin(2);
t498 = cos(t508);
t528 = pkin(2) ^ 2;
t529 = pkin(1) ^ 2;
t505 = t528 + t529;
t519 = cos(qJ(3,1));
t526 = pkin(3) ^ 2;
t556 = 0.2e1 * pkin(1);
t563 = pkin(3) * t555;
t423 = -t467 * t544 + t560;
t571 = t420 * t423;
t575 = pkin(2) * t510;
t509 = sin(pkin(7));
t576 = pkin(1) * t509;
t411 = (t417 * t519 * t563 + t420 * t526 + t429 * t505 + (pkin(3) * t417 * t498 + t429 * t575) * t556) * t527 * t467 * t429 + (t519 * t495 - t516 * t576 + pkin(3)) / (t516 * t495 + t519 * t576) * t571;
t414 = (-pkin(3) * t571 + (-t420 * pkin(3) + (-pkin(1) * t498 - t519 * pkin(2)) * t429) * t429) * t467;
t463 = mrSges(3,1) * t576 + t495 * mrSges(3,2);
t464 = t495 * mrSges(3,1) - mrSges(3,2) * t576;
t441 = t463 * t516 - t464 * t519 - Ifges(3,3);
t444 = t463 * t519 + t516 * t464;
t530 = -(m(3) * t528) - (m(2) + m(3)) * t529 - Ifges(1,3) - Ifges(2,3) - Ifges(3,3);
t537 = -t519 * mrSges(3,1) + mrSges(3,2) * t516;
t574 = m(3) * pkin(2) + mrSges(2,1);
t531 = (0.2e1 * t417 * t423 * t444 - (t537 * t555 + (-(-t537 + t574) * t510 + (mrSges(3,1) * t516 + t519 * mrSges(3,2) + mrSges(2,2)) * t509) * t556 + t530) * t414 - t441 * t411) * t467;
t579 = -t531 / 0.2e1;
t512 = legFrame(2,2);
t553 = qJ(1,2) + pkin(7);
t488 = t512 + t553;
t482 = qJ(3,2) + t488;
t489 = -t512 + t553;
t483 = qJ(3,2) + t489;
t461 = cos(t483) + cos(t482);
t501 = qJ(1,2) + t512;
t502 = qJ(1,2) - t512;
t437 = -t461 * pkin(3) + (-cos(t489) - cos(t488)) * pkin(2) + (-cos(t502) - cos(t501)) * pkin(1);
t507 = pkin(7) + qJ(3,2);
t515 = sin(qJ(3,2));
t466 = 0.1e1 / (pkin(1) * sin(t507) + t515 * pkin(2));
t431 = t437 * t466 * t550;
t458 = -sin(t482) + sin(t483);
t434 = t458 * pkin(3) + (-sin(t488) + sin(t489)) * pkin(2) + (-sin(t501) + sin(t502)) * pkin(1);
t493 = sin(qJ(1,2) + t507);
t455 = pkin(3) * t493 + pkin(2) * sin(t553) + sin(qJ(1,2)) * pkin(1);
t452 = t455 * t466 * t564;
t558 = (t458 * t567 + t461 * t566) * t466;
t569 = t493 * t521;
t416 = t431 / 0.2e1 + t452 / 0.2e1 + (t434 * t543 - t569) * t466 + t558;
t545 = t434 * t551;
t561 = t431 + t452;
t419 = (-t545 - t569) * t466 + t558 + t561;
t428 = -t466 * t569 + t558;
t497 = cos(t507);
t518 = cos(qJ(3,2));
t422 = -t466 * t545 + t561;
t572 = t419 * t422;
t410 = (t416 * t518 * t563 + t419 * t526 + t428 * t505 + (pkin(3) * t416 * t497 + t428 * t575) * t556) * t527 * t466 * t428 + (t495 * t518 - t515 * t576 + pkin(3)) / (t495 * t515 + t518 * t576) * t572;
t413 = (-pkin(3) * t572 + (-pkin(3) * t419 + (-pkin(1) * t497 - t518 * pkin(2)) * t428) * t428) * t466;
t440 = t463 * t515 - t464 * t518 - Ifges(3,3);
t443 = t463 * t518 + t515 * t464;
t538 = -t518 * mrSges(3,1) + mrSges(3,2) * t515;
t532 = (0.2e1 * t416 * t422 * t443 - (t538 * t555 + (-(-t538 + t574) * t510 + (mrSges(3,1) * t515 + t518 * mrSges(3,2) + mrSges(2,2)) * t509) * t556 + t530) * t413 - t440 * t410) * t466;
t578 = -t532 / 0.2e1;
t511 = legFrame(3,2);
t552 = qJ(1,3) + pkin(7);
t486 = t511 + t552;
t480 = qJ(3,3) + t486;
t487 = -t511 + t552;
t481 = qJ(3,3) + t487;
t460 = cos(t481) + cos(t480);
t499 = qJ(1,3) + t511;
t500 = qJ(1,3) - t511;
t436 = -t460 * pkin(3) + (-cos(t487) - cos(t486)) * pkin(2) + (-cos(t500) - cos(t499)) * pkin(1);
t506 = pkin(7) + qJ(3,3);
t514 = sin(qJ(3,3));
t465 = 0.1e1 / (pkin(1) * sin(t506) + t514 * pkin(2));
t430 = t436 * t465 * t550;
t457 = -sin(t480) + sin(t481);
t433 = t457 * pkin(3) + (-sin(t486) + sin(t487)) * pkin(2) + (-sin(t499) + sin(t500)) * pkin(1);
t492 = sin(qJ(1,3) + t506);
t454 = pkin(3) * t492 + pkin(2) * sin(t552) + sin(qJ(1,3)) * pkin(1);
t451 = t454 * t465 * t564;
t559 = (t457 * t567 + t460 * t566) * t465;
t570 = t492 * t521;
t415 = t430 / 0.2e1 + t451 / 0.2e1 + (t433 * t543 - t570) * t465 + t559;
t546 = t433 * t551;
t562 = t430 + t451;
t418 = (-t546 - t570) * t465 + t559 + t562;
t427 = -t465 * t570 + t559;
t496 = cos(t506);
t517 = cos(qJ(3,3));
t421 = -t465 * t546 + t562;
t573 = t418 * t421;
t409 = (t415 * t517 * t563 + t418 * t526 + t427 * t505 + (pkin(3) * t415 * t496 + t427 * t575) * t556) * t527 * t465 * t427 + (t495 * t517 - t514 * t576 + pkin(3)) / (t495 * t514 + t517 * t576) * t573;
t412 = (-pkin(3) * t573 + (-pkin(3) * t418 + (-pkin(1) * t496 - t517 * pkin(2)) * t427) * t427) * t465;
t439 = t463 * t514 - t464 * t517 - Ifges(3,3);
t442 = t463 * t517 + t514 * t464;
t539 = -t517 * mrSges(3,1) + mrSges(3,2) * t514;
t533 = (0.2e1 * t415 * t421 * t442 - (t539 * t555 + (-(-t539 + t574) * t510 + (mrSges(3,1) * t514 + t517 * mrSges(3,2) + mrSges(2,2)) * t509) * t556 + t530) * t412 - t439 * t409) * t465;
t577 = -t533 / 0.2e1;
t534 = (t429 ^ 2 * t444 - Ifges(3,3) * t411 + t441 * t414) * t467;
t535 = (t428 ^ 2 * t443 - Ifges(3,3) * t410 + t440 * t413) * t466;
t536 = (t427 ^ 2 * t442 - Ifges(3,3) * t409 + t439 * t412) * t465;
t1 = [t462 * t579 + t461 * t578 + t460 * t577 + (t436 * t536 + t437 * t535 + t438 * t534) * t565; t459 * t579 + t458 * t578 + t457 * t577 + (-t433 * t536 - t434 * t535 - t435 * t534) * t565; t494 * t531 + t493 * t532 + t492 * t533 + (t454 * t536 + t455 * t535 + t456 * t534) * t527;];
taucX  = t1;
