# TODO
- [ ] remove dataloader from application.hpp
- [x] skip pointcloud if distance too small -> need to change in the future, as this is discarding too many points.

# Deployment notes
- need to change number of threads in settings.hpp
- the algorithm keeps accumulating points at the same location, without smart removal of them, thus slowing down. the rate at first is 1Hz.