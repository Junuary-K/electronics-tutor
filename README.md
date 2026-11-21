# 전자공학 시각화 튜터

전자공학 전공 공부용 시각 자료 모음입니다. 각 파일은 서버 없이 단독으로 동작하는 HTML 한 개입니다.

## 도구

| 파일 | 과목 | 내용 |
|---|---|---|
| `mosfet.html` | 반도체소자 | 컷오프/선형/포화 판정, 출력·전달 특성, 채널 에너지(핀치오프), g_m·r_o |
| `mos-cap.html` | 반도체소자 | 축적/평탄/공핍/반전 판정, MOS C-V, 게이트–산화막 개략 |
| `schottky-band.html` | 반도체소자 | Schottky 밴드 다이어그램 (평형/역바이어스) |
| `three-phase.html` | 회로이론 | Y-Y·Y-Δ·Δ-Y·Δ-Δ, 평형/불평형, 4선식/3선식, 페이저도, 풀이 과정 |
| `circuit-editor.html` | 회로이론 | 회로를 직접 그리고 테브난/노턴 등가 구하기 |
| `thevenin.html` | 회로이론 | 값 편집형 테브난/노턴 단계 확인 |

`index.html` 이 목록 페이지입니다.

## 계산 검증

각 도구는 Python 스크립트로 독립 계산한 값과 대조했습니다.

- `verify_mosfet.py` — 4개 동작점, 경계 연속성, g_m·r_o 수치미분
- `verify_mosc.py` — MOS-C 7개 항목
- `verify_tyu91.py` — 교재 Schottky 예제 6개 값 + 불변식 5개
- `verify_3phase.py` — 4가지 결선 기준값
- `verify_thevenin.py` — 세 가지 경로로 구한 R_th 대조

실행: `python verify_mosfet.py` (종료 코드 0이면 통과)

## 범위와 한계

- MOSFET은 1차 square-law 모델입니다 (subthreshold, 속도포화, DIBL 없음).
- 3상은 전원이 평형일 때만 다룹니다. 부하는 상별 불평형까지 지원합니다.
- 테브난 도구는 독립전원만, 직병렬 축소만 지원합니다.

개인 학습용으로 만든 자료입니다.
