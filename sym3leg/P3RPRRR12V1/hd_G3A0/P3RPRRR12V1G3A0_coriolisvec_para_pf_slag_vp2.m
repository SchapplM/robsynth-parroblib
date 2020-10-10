% Calculate vector of centrifugal and coriolis load on the joints for
% P3RPRRR12V1G3A0
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
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d3,d4]';
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
% Datum: 2020-08-06 18:28
% Revision: 8f4ff0ee124033641e65b154ac60823cef59ef1f (2020-07-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taucX = P3RPRRR12V1G3A0_coriolisvec_para_pf_slag_vp2(xP, xDP, qJ, legFrame, ...
  koppelP, pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,1),zeros(3,3),zeros(3,3),zeros(3,3),zeros(6,1),zeros(3+1,1),zeros(3+1,3),zeros(3+1,6)}
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3RPRRR12V1G3A0_coriolisvec_para_pf_slag_vp2: xP has to be [3x1] (double)');
assert(isreal(xDP) && all(size(xDP) == [3 1]), ...
  'P3RPRRR12V1G3A0_coriolisvec_para_pf_slag_vp2: xDP has to be [3x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 3]), ...
  'P3RPRRR12V1G3A0_coriolisvec_para_pf_slag_vp2: qJ has to be [3x3] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'P3RPRRR12V1G3A0_coriolisvec_para_pf_slag_vp2: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [4 1]), ...
  'P3RPRRR12V1G3A0_coriolisvec_para_pf_slag_vp2: m has to be [4x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [4,3]), ...
  'P3RPRRR12V1G3A0_coriolisvec_para_pf_slag_vp2: mrSges has to be [4x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [4 6]), ...
  'P3RPRRR12V1G3A0_coriolisvec_para_pf_slag_vp2: Ifges has to be [4x6] (double)'); 
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3RPRRR12V1G3A0_coriolisvec_para_pf_slag_vp2: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3RPRRR12V1G3A0_coriolisvec_para_pf_slag_vp2: Koppelpunkt has to be [3x3] (double)');

%% Symbolic Calculation
% From coriolisvec_para_plfcoord_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2020-08-06 18:28:08
% EndTime: 2020-08-06 18:28:10
% DurationCPUTime: 1.69s
% Computational Cost: add. (6516->222), mult. (9714->377), div. (2883->7), fcn. (8442->18), ass. (0->149)
t502 = legFrame(3,2);
t489 = cos(t502);
t511 = cos(qJ(3,3));
t497 = t511 ^ 2;
t505 = sin(qJ(3,3));
t506 = sin(qJ(1,3));
t585 = (-pkin(5) - pkin(6));
t492 = pkin(1) - t585;
t512 = cos(qJ(1,3));
t535 = qJ(2,3) * t506 + t492 * t512;
t486 = sin(t502);
t570 = t486 * t511;
t438 = t535 * t489 * t505 + qJ(2,3) * t570 + (t505 * t570 + (-t497 + 0.1e1) * t489 * t506) * pkin(3);
t567 = t489 * t511;
t441 = (pkin(3) * t567 - t535 * t486) * t505 + t506 * pkin(3) * (t511 - 0.1e1) * (t511 + 0.1e1) * t486 + qJ(2,3) * t567;
t475 = pkin(3) * t505 + qJ(2,3);
t465 = t475 * t512 - t492 * t506;
t472 = 0.1e1 / t475;
t494 = 0.1e1 / t505;
t517 = xDP(3);
t518 = xDP(2);
t519 = xDP(1);
t594 = 0.2e1 * (t465 * t517 + (t438 * t519 + t441 * t518) * t494) * t472;
t503 = legFrame(2,2);
t490 = cos(t503);
t513 = cos(qJ(3,2));
t498 = t513 ^ 2;
t507 = sin(qJ(3,2));
t508 = sin(qJ(1,2));
t514 = cos(qJ(1,2));
t536 = qJ(2,2) * t508 + t492 * t514;
t487 = sin(t503);
t569 = t487 * t513;
t439 = t536 * t490 * t507 + qJ(2,2) * t569 + (t507 * t569 + (-t498 + 0.1e1) * t490 * t508) * pkin(3);
t566 = t490 * t513;
t442 = (pkin(3) * t566 - t536 * t487) * t507 + t508 * pkin(3) * (t513 - 0.1e1) * (t513 + 0.1e1) * t487 + qJ(2,2) * t566;
t476 = pkin(3) * t507 + qJ(2,2);
t466 = t476 * t514 - t492 * t508;
t473 = 0.1e1 / t476;
t495 = 0.1e1 / t507;
t593 = 0.2e1 * (t466 * t517 + (t439 * t519 + t442 * t518) * t495) * t473;
t504 = legFrame(1,2);
t491 = cos(t504);
t515 = cos(qJ(3,1));
t499 = t515 ^ 2;
t509 = sin(qJ(3,1));
t510 = sin(qJ(1,1));
t516 = cos(qJ(1,1));
t537 = qJ(2,1) * t510 + t492 * t516;
t488 = sin(t504);
t568 = t488 * t515;
t440 = t537 * t491 * t509 + qJ(2,1) * t568 + (t509 * t568 + (-t499 + 0.1e1) * t491 * t510) * pkin(3);
t565 = t491 * t515;
t443 = (pkin(3) * t565 - t537 * t488) * t509 + t510 * pkin(3) * (t515 - 0.1e1) * (t515 + 0.1e1) * t488 + qJ(2,1) * t565;
t477 = pkin(3) * t509 + qJ(2,1);
t467 = t477 * t516 - t492 * t510;
t474 = 0.1e1 / t477;
t496 = 0.1e1 / t509;
t592 = 0.2e1 * (t467 * t517 + (t440 * t519 + t443 * t518) * t496) * t474;
t447 = (-t506 * t517 + (-t486 * t518 + t489 * t519) * t512) * t472;
t528 = 0.1e1 / pkin(3);
t459 = (-t486 * t519 - t489 * t518) * t528 * t494;
t522 = qJ(2,3) ^ 2;
t527 = pkin(3) ^ 2;
t526 = pkin(5) ^ 2;
t529 = pkin(1) ^ 2;
t588 = -2 * pkin(5);
t534 = (2 * t585 * pkin(1)) - t526 - t527 - t529 + ((t588 - pkin(6)) * pkin(6));
t584 = pkin(3) * t459;
t555 = t511 * t584;
t577 = qJ(2,3) * t505;
t556 = -0.2e1 * t577;
t576 = t447 * t492;
t417 = (((t511 * t576 - t584) * t505 - qJ(2,3) * t459) * t494 * t584 + t576 * t594 + (t492 * t555 + (pkin(3) * t556 + t497 * t527 - t522 + t534) * t447) * t447) * t472;
t420 = (t594 + 0.2e1 * t555 - t576) * t447 * t472;
t456 = t459 ^ 2;
t520 = pkin(1) + pkin(5);
t481 = mrSges(3,2) * t520 - Ifges(3,6);
t482 = mrSges(3,1) * t520 - Ifges(3,5);
t462 = -t481 * t505 + t482 * t511;
t580 = (mrSges(3,3) - mrSges(2,2));
t471 = m(2) * pkin(1) + m(3) * t520 + t580;
t521 = m(2) + m(3);
t478 = qJ(2,3) * t521 + mrSges(2,3);
t501 = Ifges(3,1) - Ifges(3,2);
t530 = -m(3) * t526 + (mrSges(3,3) * t588) - t521 * t529 - Ifges(2,1) - Ifges(3,2) - Ifges(1,3);
t531 = (mrSges(3,1) * qJ(2,3) - t501 * t505) * t511;
t543 = -mrSges(3,1) * t505 - mrSges(3,2) * t511;
t549 = t456 * t494 * t511;
t552 = mrSges(3,2) * t577;
t581 = (m(3) * pkin(5) + t580) * pkin(1);
t587 = -0.2e1 * mrSges(2,3);
t564 = (mrSges(3,1) * t556 + qJ(2,3) * t587 - t501 * t497 - t521 * t522 + t530) * t420 + t471 * t417 + t462 * t549 + 0.2e1 * (-(mrSges(3,2) * qJ(2,3) - Ifges(3,4) * t505) * t511 - t581) * t420 + (t481 * t511 + t482 * t505) * t456 + ((t478 - t543) * t594 + (-0.2e1 * t552 + (-0.4e1 * t497 + 0.2e1) * Ifges(3,4) + 0.2e1 * t531) * t459) * t447;
t591 = t512 * t564;
t448 = (-t508 * t517 + (-t487 * t518 + t490 * t519) * t514) * t473;
t460 = (-t487 * t519 - t490 * t518) * t528 * t495;
t523 = qJ(2,2) ^ 2;
t583 = pkin(3) * t460;
t554 = t513 * t583;
t578 = qJ(2,2) * t507;
t557 = -0.2e1 * t578;
t575 = t448 * t492;
t418 = (((t513 * t575 - t583) * t507 - qJ(2,2) * t460) * t495 * t583 + t575 * t593 + (t492 * t554 + (pkin(3) * t557 + t498 * t527 - t523 + t534) * t448) * t448) * t473;
t421 = (t593 + 0.2e1 * t554 - t575) * t448 * t473;
t457 = t460 ^ 2;
t463 = -t481 * t507 + t482 * t513;
t479 = qJ(2,2) * t521 + mrSges(2,3);
t532 = (mrSges(3,1) * qJ(2,2) - t501 * t507) * t513;
t542 = -mrSges(3,1) * t507 - mrSges(3,2) * t513;
t548 = t457 * t495 * t513;
t551 = mrSges(3,2) * t578;
t563 = (mrSges(3,1) * t557 + qJ(2,2) * t587 - t501 * t498 - t521 * t523 + t530) * t421 + t471 * t418 + t463 * t548 + 0.2e1 * (-(mrSges(3,2) * qJ(2,2) - Ifges(3,4) * t507) * t513 - t581) * t421 + (t481 * t513 + t482 * t507) * t457 + ((t479 - t542) * t593 + (-0.2e1 * t551 + (-0.4e1 * t498 + 0.2e1) * Ifges(3,4) + 0.2e1 * t532) * t460) * t448;
t590 = t514 * t563;
t449 = (-t510 * t517 + (-t488 * t518 + t491 * t519) * t516) * t474;
t461 = (-t488 * t519 - t491 * t518) * t528 * t496;
t524 = qJ(2,1) ^ 2;
t582 = pkin(3) * t461;
t553 = t515 * t582;
t579 = qJ(2,1) * t509;
t558 = -0.2e1 * t579;
t574 = t449 * t492;
t419 = (((t515 * t574 - t582) * t509 - qJ(2,1) * t461) * t496 * t582 + t574 * t592 + (t492 * t553 + (pkin(3) * t558 + t499 * t527 - t524 + t534) * t449) * t449) * t474;
t422 = (t592 + 0.2e1 * t553 - t574) * t449 * t474;
t458 = t461 ^ 2;
t464 = -t481 * t509 + t482 * t515;
t480 = qJ(2,1) * t521 + mrSges(2,3);
t533 = (mrSges(3,1) * qJ(2,1) - t501 * t509) * t515;
t541 = -mrSges(3,1) * t509 - mrSges(3,2) * t515;
t547 = t458 * t496 * t515;
t550 = mrSges(3,2) * t579;
t562 = (mrSges(3,1) * t558 + qJ(2,1) * t587 - t501 * t499 - t521 * t524 + t530) * t422 + t471 * t419 + t464 * t547 + 0.2e1 * (-(mrSges(3,2) * qJ(2,1) - Ifges(3,4) * t509) * t515 - t581) * t422 + (t481 * t515 + t482 * t509) * t458 + ((t480 - t541) * t592 + (-0.2e1 * t550 + (-0.4e1 * t499 + 0.2e1) * Ifges(3,4) + 0.2e1 * t533) * t461) * t449;
t589 = t516 * t562;
t586 = -0.2e1 * Ifges(3,4);
t444 = t447 ^ 2;
t468 = -mrSges(3,1) * t511 + mrSges(3,2) * t505;
t561 = -t417 * t521 + t420 * t471 + t468 * t549 - t444 * t478 + t543 * (t444 + t456);
t445 = t448 ^ 2;
t469 = -mrSges(3,1) * t513 + mrSges(3,2) * t507;
t560 = -t418 * t521 + t421 * t471 + t469 * t548 - t445 * t479 + t542 * (t445 + t457);
t446 = t449 ^ 2;
t470 = -mrSges(3,1) * t515 + mrSges(3,2) * t509;
t559 = -t419 * t521 + t422 * t471 + t470 * t547 - t446 * t480 + t541 * (t446 + t458);
t546 = t561 * t494;
t545 = t560 * t495;
t544 = t559 * t496;
t540 = (Ifges(3,3) * t549 - t417 * t468 - t420 * t462 + t444 * (t497 * t586 + Ifges(3,4) + t531 - t552)) * t494;
t539 = (Ifges(3,3) * t548 - t418 * t469 - t421 * t463 + t445 * (t498 * t586 + Ifges(3,4) + t532 - t551)) * t495;
t538 = (Ifges(3,3) * t547 - t419 * t470 - t422 * t464 + t446 * (t499 * t586 + Ifges(3,4) + t533 - t550)) * t496;
t1 = [(t440 * t544 + t491 * t589) * t474 + (t439 * t545 + t490 * t590) * t473 + (t438 * t546 + t489 * t591) * t472 + (t486 * t540 + t487 * t539 + t488 * t538) * t528; (t443 * t544 - t488 * t589) * t474 + (t442 * t545 - t487 * t590) * t473 + (t441 * t546 - t486 * t591) * t472 + (t489 * t540 + t490 * t539 + t491 * t538) * t528; (t559 * t467 - t562 * t510) * t474 + (t560 * t466 - t563 * t508) * t473 + (t561 * t465 - t564 * t506) * t472;];
taucX  = t1;
